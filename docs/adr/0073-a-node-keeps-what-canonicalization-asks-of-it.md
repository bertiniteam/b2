# ADR-0073: A node keeps what canonicalization asks of it; differentiation visits each node once

**Status:** Accepted
**Date:** 2026-10-11

## Context

Building a sum or a product puts its operands in canonical order (ADR-0059), and to do that it asks each operand three things: its variables, its multidegree, and, when two multidegrees tie, its canonical encoding.  Each answer was a walk of the operand.  #417/#418 made the multidegree walk visit each node once per call, so a question cost the operand's DAG rather than its expansion, and recorded that a cache of results between calls had measured 29% slower, because filling it cost the same walk.

That was still a walk per construction.  Differentiation builds as many nodes as the expression has, each over operands it has just built, so building a derivative cost the square of its size.  Two holes made it worse:

- `Differentiate` had no memo at all.  A subexpression reached by many paths was differentiated once per path: the sub-minors of a determinant built by cofactor expansion, k! times at the top.
- A unary operator's multidegree is computed from `Degree(v)` of its operand, and `Degree(v)` walked the whole expansion.  So every construction above a negation reopened the walk #418 had closed.

Deflation appends Jacobian minors, so these costs landed on every singular point a decomposition sharpens: on the Barth sextic, 121 s of first evaluations for 50 deflated systems, against 0.1 s of evaluating them (#505).

## Decision

- **`Differentiate` is non-virtual and keeps each derivative for the outermost call.**  A node type's rule is `DifferentiateImpl`, which differentiates its operands by calling `Differentiate`.  A `node::DerivativeScope` extends the memo across a loop; both Jacobian builds open one, so functions that share subexpressions differentiate each of them once among them.
- **Each node keeps a degree record, made once from its operands' records.**  It holds the node's variables (what `GatherVariables` returns), its degree and multidegree entry in each variable it mentions (its variables and its differentials' variables), its total degree, and one degree and one multidegree entry for every variable it does not mention.  `Degree`, `MultiDegree` and `GatherVariables` read it.  It is made without a walk, so #418's measurement does not apply to it.
- **The answer for an unmentioned variable is computed, not assumed.**  It is not always 0: a power with a negative constant exponent reports -1 in every variable.  The degree rules compare a variable only by identity with the variables a subtree mentions, so every unmentioned variable gets the same answer, and the record computes it once by asking the rules about a variable that is in no expression (a default-constructed `Variable`, never interned, so no user's name can collide with it).
- **A node keeps its canonical encoding once computed**, as the tie-break key.  It is the same text ADR-0059 chose; only where it is kept changes.
- **Records and keys are published once, by atomic compare-and-swap.**  Systems compile lazily on worker threads, so two threads may make the same record; the first one published wins, and a reference into it lives as long as the node.

## Consequences

- A node type's `DegreeImpl` and `MultiDegreeImpl` must never ask `Degree` or `MultiDegree` of their own node: the record is being made from them, and the request would recurse until the stack overflows.  Ask `DegreeImpl` of `this` instead, as the unary, power and named-expression rules do.
- Do not replace the unmentioned answer with 0, and do not key the sentinel variable by a name.
- Do not go back to walking an operand to put a sum or product in order, or to a memo that lasts one call: a construction's question should cost O(operands).
- Do not make `Differentiate` or `Degree(v)` virtual again; override the `Impl`.
- The canonical order is unchanged, so no digest moves and no encoding version is bumped.  The canonical encoding itself still costs a walk of the node it encodes, once per node, because its back-references number the whole tree; it is asked only on a multidegree tie.
- Memory: every node that has been an operand of a canonicalized sum or product holds a record, proportional to the variables it mentions, and possibly its encoding.
