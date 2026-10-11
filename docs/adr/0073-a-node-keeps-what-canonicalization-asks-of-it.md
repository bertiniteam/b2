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
- **The tie-break key is a content key, a Merkle hash of the canonical encoding** (this replaces ADR-0059's key, the operand's whole canonical encoding).  A node's key is the SHA-256 of its shallow canonical encoding: its own kind, payload and operands in the ADR-0042 format, each operand written as `@` and that operand's key.  It is made from the operands' keys, once per node, and kept.  The whole encoding could not be: its back-references number the whole tree, so it is not made from the operands' encodings, and computing it once per node still walked each node's DAG -- O(DAG^2) to build an expression, 1.4 s for the derivative of a 7 x 7 determinant against 0.04 s with the key.  The key is still derived from the encoding, as ADR-0059 asks of any ordering criterion, exact, and the same on every platform.
- **Degree records and content keys are published once, by atomic compare-and-swap.**  Systems compile lazily on worker threads, so two threads may make the same record; the first one published wins, and a reference into it lives as long as the node.

## Consequences

- A node type's `DegreeImpl` and `MultiDegreeImpl` must never ask `Degree` or `MultiDegree` of their own node: the record is being made from them, and the request would recurse until the stack overflows.  Ask `DegreeImpl` of `this` instead, as the unary, power and named-expression rules do.
- Do not replace the unmentioned answer with 0, and do not key the sentinel variable by a name.
- Do not go back to walking an operand to put a sum or product in order, or to a memo that lasts one call: a construction's question should cost O(operands).
- Do not make `Differentiate` or `Degree(v)` virtual again; override the `Impl`.
- Where operands tie on multidegree, their canonical order changes, so the system encoding version is `b2sysenc/3` and every system digest moved in 4.0.1: records written by 4.0.0 are asked again, not recalled.  The golden fixture has a recipe whose operands tie (`tied_terms`), so a later change to the tie-break moves the registry's keyspace hash; before it, no recipe had a tie and the registry could not see one.
- Do not order operands by anything that walks them, the full encoding included.  A new tie-break must be made from the operands' keys, or it is a walk again.
- A collision of two different operands' keys would leave them in input order: a valid expression, not a canonical one.  At SHA-256 that is not a practical concern.
- Memory: every node that has been an operand of a canonicalized sum or product holds a record, proportional to the variables it mentions, and, if it was ever in a tie, its 64-character content key.
