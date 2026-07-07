# ADR-0028: Named-node taxonomy — NameHolder, three named kinds, NamedExpression replaces Handle, Jacobian deleted

**Status:** Accepted (supersedes the "keep `Handle` for `Jacobian`" decision in ADR-0027)
**Date:** 2026-06-22

> **Implementation status (2026-06-22, on `feature/eliminate_duplicate_nodes` @ `c26ea382`):**
> *Done:* `NameHolder` rename; `NamedExpression` + `Named(expr,name)`; `Find<T>` (`GatherVariables`
> = `Find<Variable>`); **`node::Jacobian` deleted** (`b3573117`); **subfunctions → `NamedExpression`**
> (`c26ea382`) — storage deleted, `Describe` discovers via `Find<NamedExpression>`, and `FunctionParser`
> no longer wraps results in a `Function` (clean printing).
> *Correction to this ADR:* `differentiate()` already returns a **`Differential` tree, not a dict** — the
> "`{variable: derivative}` dict" below was a misread of the code; **no dict change was made or is needed.**
> *Remaining = **Track F** (`z_notes/2026-06-22_roadmap_FE.md`):* **F1** rewrite the classic parser to
> eager-bind the declare-then-define kinds (`function`/`constant`/`parameter`) + convert their System
> storage (`explicit_parameters_`, `constant_subfunctions_`) to `NamedExpression` — then **F2** delete the
> `Function` class and fold `Handle` into `NamedExpression`. Detailed steps in
> `z_notes/2026-06-22_finish_function_handle_removal.md`. (An F1 draft was reverted to keep the branch green.)

## Context

Eliminating `node::Function` (ADR-0027) raised the question of whether its base, `Handle`, is a
wasted layer. For *top-level functions* it was: the forwarding is redundant (a bare expression node
does everything itself), the name was vestigial, and the mutability (`SetRoot`) was a liability the
hash-consing/immutability design wants gone. Those were removed — parsed functions are now eager-bound
bare `Nd`.

But for *named subexpressions* (`a = x^2 + y^2`, then `a^2 + a + 1`), the name is genuinely **used**
(you want to print `a`, and see `a = x^2+y^2` below) and the node is genuinely **embedded** (other
expressions reference `a`). So `Handle` is not wasted there — it is the *named expression* concept.
What was wrong with it was the **mutability** and the **special storage**, not the idea.

## Decision

### 1. `NameHolder` mixin

Rename the current `Named` name-storage mixin to **`NameHolder`** (an `IsNamed` concept if/when we
adopt C++20). It is the shared "this node has a name" capability.

### 2. Three named *kinds* — distinct concrete types, each with its own `Find`

- **`NamedSymbol`** — π, e: names **baked in by the core**, not user-chosen.
- **`Variable`** — user-named leaves.
- **`NamedExpression`** — user-named wrappers around an expression (the surviving, honest form of
  `Handle`).

All three are `NameHolder`s; nothing else is. They are *separate concrete types*, so discovery
discriminates: `Find(Variable)`, `Find(NamedSymbol)`, `Find(NamedExpression)` each return only their
own kind. `Find(NamedExpression)` never returns variables or π/e.

### 3. `NamedExpression` (replaces `Function`/`Handle`)

- **Immutable / eager**: `Named(expr, "a")` takes the expression **at construction**. No `SetRoot`,
  no declare-then-fill. The parser builds the expression first, then names it.
- **Hash-consed** by (expression, name): `Named(e,"a")` is distinct from bare `e` and from
  `Named(e,"b")`. It evaluates to its expression (one copy in the SLP); the inner expression is still
  shared/computed-once by hash-consing, so the wrapper costs ~nothing.
- **Prints as its name.** The expansion is revealed elsewhere (`Describe`), which lists `a = <expr>`
  for each discovered `NamedExpression`. Nested names (`b = a + 1`) fall out of recursive `Find`.
- **No special storage.** Named subexpressions just live in the trees and are **discovered**.
  `block.subfunctions_` / `constant_subfunctions_` / `AddSubFunction` / `AddConstant` /
  `ConstantSubfunctions` / `NumConstants` are **deleted**.

### 4. "Function" is a *role*, not a kind

A function is a bare `Nd` output the System evaluates — optionally a `NamedExpression`. Python
functions are unnamed by default (`add_function(x*y)`); the Bertini-1 parser, which must name them,
wraps each in a `NamedExpression`. So top-level function names return — as `NamedExpression` display,
not as the dead `name` string parameter that ADR-0027's cut removed (this is the better direction:
`Describe` can show the user's `f`, not `f_0`).

### 5. `Find` (rename of `Gather`)

Rename `Gather` → `Find` (sympy's `find`). `Find(T)` returns the distinct nodes of kind `T` in a
subtree (deduped by identity, name-ordered). `GatherVariables` becomes `Find(Variable)`.

### 6. `node::Jacobian` is deleted (and `Handle` with it)

`node::Jacobian` is a legacy `Handle` subclass carrying a mutable `current_diff_variable_`, evaluated
by the old tree-walk `EvalJ`. **The SLP never compiled it** (`SLPCompiler::Visit(node::Jacobian)`
throws "unimplemented"); production derivatives are per-variable bare `Differentiate(var)` nodes. Its
only live use is the Python all-variables `differentiate()`. So it is **deleted**, not collapsed into
anything. The all-variables `differentiate()` returns a **dict `{variable: derivative}`** (sympy
flavor); the single-variable `differentiate(var)` already returns a plain `Node`. With `Jacobian`
gone, `Handle`'s last justification is gone — `Handle` is deleted, replaced wholly by
`NamedExpression`.

## Consequences

- The symbol layer becomes a clean trichotomy of named kinds over one `NameHolder` mixin; "function"
  stops being a node type and becomes an output role.
- Build order (each a green increment where possible): `NameHolder` rename → introduce immutable
  `NamedExpression` + `Named(expr, name)` (C++/Python) → `Find` (rename `Gather`, add per-kind
  finders) → repoint the parser to make `NamedExpression`s + delete subfunction storage + `Describe`
  discovers → delete `node::Jacobian` + `Handle`, `differentiate()` → dict.
- Supersedes ADR-0027 §"keep `Handle` (base of `node::Jacobian`)"; the rest of ADR-0027 (output-set
  compilation, Program/Memory split, freeze-set partition) stands.

References: ADR-0027; this branch's `6b837211` (compiler bare roots), `8ac65cd0` (parser eager-bind),
`fd93dbff` (functions stored bare). The `constant` keyword is obsoleted by the freeze-set partition
(a constant is auto-frozen), so it degrades to a plain `NamedExpression`.
