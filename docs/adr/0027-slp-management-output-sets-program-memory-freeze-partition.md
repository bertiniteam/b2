# ADR-0027: SLPs belong to output-sets, not nodes; Program/Memory split with a freeze-set tape partition

**Status:** Accepted (design; implementation staged)
**Date:** 2026-06-22

> **Implementation status (2026-06-22, on `feature/eliminate_duplicate_nodes` @ `c26ea382`):**
> *Built:* SLP is the sole evaluator (`942f7f18`); eval-without-System **adapter** (`a052de37`,
> a stepping stone); freeze-set tape partition (the per-step constant-folding prologue);
> output-sets compiled from explicit bare roots. *Revised by ADR-0028:* `node::Function`'s name
> role is carried by **`NamedExpression`**, not a parallel `vector<string>` as written below;
> and `Handle` does **not** survive (no remaining subclass once `Jacobian` was deleted).
> **E1 DONE (2026-06-22, `f9ea0258`..`a5a17d3b`):** the SLP is split into an immutable, shareable
> `SLPProgram` (instructions + exact node-free `ConstantRecipe`s + layout) and a per-thread
> `SLPMemory` (registers + precision + freshness flags); the per-thread eval path writes **no**
> shared node state (variables/path-value/precision are System-owned buffers + the SLP's own
> Memory, all node-free); and `Clone(System)` is now a **Memory-isolating shallow copy** — it
> shares the immutable DAG + compiled Program and copies only the eval Memory (operand-holding
> blocks deep-copy their nested operand Systems the same way), deleting the **#246** serialize/
> recompile deep copy. Caveat: node value/precision *storage* still exists (E5 deletes it); a
> temporary `ConstantRecipe::Snapshot` kind handles fixed-variables-as-constants until then.
> *Remaining **Track E** (see `z_notes/2026-06-22_roadmap_FE.md`):* **E2** retire the adapter
> (real `Compile([outputs],…)` + on-node Program memo) → **E5** delete node-level eval [DEFERRED];
> **E3** currying/freeze-set API → **E6** tutorial; **E4** intern identical Programs — **DONE**
> (ADR-0042: `SLPProgram::ContentHash`/`SameContent` + `InternProgram` weak table, wired into
> `SLPCompiler::Compile` and facade deserialization).

## Context

The `function_tree` refactor (hash-consing in PR #25, then SLP-as-sole-evaluator and
eval-without-System on `feature/slp_only_eval`) keeps converging on one unanswered question:
**what owns a compiled SLP, and which nodes are evaluation entry points?**

The history of `node::Function` is the history of that confusion. `Function` began as a *name*
for an equation. It then acquired a second job: marking *which nodes are top-level outputs* so the
SLP compiler knows where to write results (`SLPCompiler` keys output slots on `Function`
pointers, `locations_top_level_functions_and_derivatives_`). But SLPs are compiled from
**Systems**, and a System already *is* its list of top-level nodes. So the entry-point role was
always redundant with "the caller knows its own outputs," and `Function` is, stripped of its name
role, useless.

Two failure modes haunt any fix:

- **#nodes-many SLPs.** If a Program is compiled per node (e.g. by routing `Node::Eval` through a
  throwaway per-call SLP, or by caching a Program on *every* node), then `a+b+c` yields five SLPs
  (`a`, `b`, `c`, `a+b`, `a+b+c`). That is absurd: only the root is ever evaluated.
- **Recompile-per-call.** Compiling a fresh SLP on every `Eval` is a leak, not an engine.

Separately, the project wants **optimal evaluation performance**, and the standing observation is
that *constant sub-expressions* (the γ, coefficient products, `2·pi`) need recomputing only when
**precision** changes, never when the **point** changes — and users must not have to find them by
hand. The old `Reset`/`SetVariableValues` gestured at this but never delivered, because the
tree-walk's down-only invalidation could not reliably assert "this subtree is constant, leave it."

## Decision

### 1. The unit that owns an SLP is an *output-set*, never a node.

Recast the one compile operation to take outputs explicitly:

    SLPCompiler::Compile(outputs: vector<Nd>, variable_order) -> Program
        // one output slot per root Nd, in order; all shared internals CSE'd exactly once

**Entry points are the explicit `outputs` argument — full stop.** A System passes its functions
(and their derivatives); a bare `f.eval` passes `[f]`. Output slots are keyed by the **root `Nd`**,
not by a `Function` pointer.

This kills both failure modes structurally. Handed root `a+b+c`, the compiler walks the DAG **once**,
gives `a`, `b`, `c` internal register slots, and emits **one** Program with one output — `a+b` is an
*instruction*, not a sub-Program. The number of Programs that exist equals the number of output-sets
actually compiled, never the node count.

`node::Function` is **eliminated.** Its name role moves to parallel `vector<string>` metadata on the
System; its entry-point role becomes the `outputs` argument. (`Handle` survives as the base of the
symbolic `node::Jacobian`.)

### 2. A bare expression memoizes its one-output Program on the node.

`f.eval(...)` lazily compiles the Program for `([f], canonical-by-name order)` and caches it on the
**immutable, hash-consed** node — exactly like the memoized `Hash()`. Hash-consing makes the cache
auto-shared (identical expressions are the same node) and never stale (nodes never change). This does
**not** reintroduce #nodes-many SLPs: only a node actually evaluated *as a root* ever gets a cached
Program, and caching it does not compile its children (one root compile = one whole-subtree Program).

### 3. Program / Memory split.

- **Program** — instructions + the constant *recipe* (precision-free). Immutable, shareable
  **read-only across threads.**
- **Memory** — the register file + the working precision + the current point. **Per-thread.**

The honest prize is **not** "`System::Clone()` vanishes" — path tracking still does per-thread
mutable work. The prize is that the **deep clone of the node DAG and block operands (#246)** is
deleted: threads share the immutable DAG + Program and allocate a Memory (cheap) instead of
deep-copying the whole tree. Whether a featherweight per-thread System *view* remains is an
implementation detail, not a guarantee. The enabler already exists: because the SLP is the sole
evaluator, the per-thread point already lives in the SLP's own register file (`SetVariableValues`
copies it in), and `Variable::current_value_` is now vestigial during tracking (production no longer
node-evaluates). What remains is to move *precision* into Memory and stop writing vestigial node
state.

### 4. Freeze-set tape partition — generalized from day one.

The compiler partitions the instruction tape by **input-dependency**, against a **freeze set** of
inputs:

- the **frozen segment** — instructions reachable only from frozen inputs; recomputed only when a
  frozen input (or precision) changes;
- the **live segment** — instructions with a live (non-frozen) ancestor; recomputed when the live
  inputs change.

This is **one** parameterized mechanism, not two:

- **constant-folding = freeze set `{}`** (only literals are frozen). The constant prologue is
  recomputed only on **precision** change — the every-step tracker win, since γ/coefficients are
  otherwise recomputed every step for nothing.
- **currying `f(x=…)` over varying `y` = freeze set `{x}`** — the x-only sub-tape becomes a second
  frozen prologue, recomputed only when `x` (or precision) changes.

"Is this register constant?" is `subtree contains no Variable`, a memoizable boolean on the
hash-consed DAG (like `Degree()`), computed by the **compiler** — **zero user burden**. Memory gates
each frozen segment with a stamp (precision for the constant prologue; precision + frozen values for
a curried segment). Making the freeze set a first-class compiler/Program parameter from the start
means constant-folding is its default instance and currying needs no later redesign.

### 5. No output-freshness flags.

Per-output "already computed f this step" flags are **dropped.** They exist only to skip a redundant
*same-point* re-eval — which assumes a careless caller and reintroduces the stale-cache footgun the
whole refactor is killing. In tracking they barely fire anyway (every Newton iterate is a new point;
predictor and corrector evaluate at different points). Evaluation is **recompute-on-demand**; the
frozen-segment skip already makes "recompute" cheap. Trust the caller. Freshness is *not* the
mechanism behind constant-folding (that is the tape partition) and *not* behind currying (the same
partition with a non-empty freeze set), so nothing needs per-output or per-node freshness state.

## Consequences

- One coherent refactor subsumes three threads that looked separate: eval-engine unification,
  `Function`-elimination, and the Program/Memory split. `Compile(outputs, vars)` is the keystone.
- Memory holds exactly: registers + working precision + current point + a stamp per frozen segment.
  No per-output, no per-node freshness.
- The adapter `EvalExpression` shipped on this branch (one-function throwaway System per call) is a
  **stepping stone**; it is replaced by `Compile([f], canonical)` + the on-node Program memo, behind
  the already-frozen `f.eval` / `EvalExpression` API.
- Staged, each a green increment where possible:
  1. **Freeze-set tape partition inside the current SLP** (independent of the rest): mark constant
     registers (memoized "no Variable descendant"), emit a constant prologue, skip it on
     point-only changes; structure the partition around a freeze set with `{}` as the default.
     Delivers the per-step perf win immediately; verify against `benchmark/baseline_*.csv`.
  2. **`Compile(outputs: vector<Nd>, vars)`** + re-key output slots on `Nd`. Interdependent with
     `Function`-elimination (no safe mid-commits): one focused push to green.
  3. Replace the bare-eval adapter with `Compile([f], canonical)` + on-node Program memo.
  4. **Program/Memory split**; move precision into Memory; stop writing vestigial node values;
     delete #246's deep operand clone. Parallel crossed-paths tracking shares one Program.
  5. Widen the freeze set through the public API when currying lands (no redesign).
- Deferred: output-slicing (evaluate only the requested outputs' instructions) is an *orthogonal*
  tape optimization, separate from input-freezing.

References: this branch's `942f7f18` (SLP sole evaluator), `a052de37` (eval-without-System adapter),
`1bb8e17e` (blend coefficients via System, last production tree-eval removed); ADR-0011
(immutability), issues #246, #251; the plan note `z_notes/2026-06-21_eval_engine_refactor_plan.md`.
