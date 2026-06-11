# ADR-0010: function_tree uses single inheritance; cross-cutting traits are capability classes

**Status:** Accepted
**Date:** 2026-06-11 (PR #6)

## Context

Every class in `function_tree/` inherited with `public virtual`, forced by exactly two
diamond sources:

1. `EnableSharedFromThisVirtual<T> : public virtual Node` — a mixin giving derived
   classes a typed `shared_from_this`. Mixed into ~24 concrete classes, it created a
   diamond with each class's real base chain (its own comment said it existed to
   "solve the diamond problem" — the one it created). It was barely load-bearing:
   every SLP call site already did its own `dynamic_pointer_cast` on
   `Node::shared_from_this()`; the single typed caller was `Variable::Differentiate`.
2. `Pi`/`E : public virtual Number, public virtual NamedSymbol` — the only genuine
   semantic diamond (both bases derive from `Symbol`).

Costs of the virtual inheritance: dynamic base-pointer adjustments in the eval hot
paths, awkward `dynamic_pointer_cast` where a static cast should do, dual-path
serialization for Pi/E (`base_object<Number>` *and* `base_object<NamedSymbol>` into a
shared virtual base), constructor rules that let derived classes initialize
grandparent bases (six trig ctors did), and a Python-visible hierarchy that could not
be expressed in single-inheritance binding frameworks (nanobind, a likely future
target — see the migration assessment notes).

Which base do Pi/E keep? Differentiation's constant detection does
`dynamic_pointer_cast<Number>` (`SumOperator::Differentiate`), and Pi/E flow through
it — so they must remain `Number`s. Nothing in core ever casts to `NamedSymbol`; its
only contribution is a `name_` string and a one-line `print`.

## Decision

- **Delete the mixin.** `Node` keeps plain `std::enable_shared_from_this<Node>`; the
  one typed caller casts explicitly (`static_pointer_cast`, valid under non-virtual
  inheritance).
- **Cross-cutting traits become capability classes that do NOT derive from Node.**
  `Named` (symbol.hpp) carries `name_`/`name()`; `Pi : public Number, public Named`.
  A capability class adds no second path to the Node root, needs no vtable (protected
  non-virtual destructor), and costs nothing in the bound hierarchy.
- **De-virtualize the whole subsystem**: all `public virtual` base specifiers in
  function_tree became plain `public`, including `Node : VisitableBase<>` (its sole
  inheritor). Derived constructors initialize only their direct base (trig ctors call
  `TrigOperator(N)`, which forwards).
- Python bindings: Pi/E declare `bases<Number>` and apply `NamedSymbolVisitor`
  directly (they no longer descend from `AbstractNamedSymbol`), so `pi.name` is
  unchanged from Python.

## Consequences

**Positive**
- Static base adjustments; plain `shared_from_this`; single-path serialization;
  every bound class is straightforwardly single-inheritance (nanobind-ready).
- The hierarchy is honest: one IS-A chain per class, capabilities by composition.

**Negative / to watch**
- Object layout and Pi/E archive layout changed: anything linking the old
  `libbertini2` needs a rebuild (automatic in-repo; no archives persist across
  versions).
- New classes must NOT reintroduce a second Node-derived base. If a node needs a
  cross-cutting trait, write a capability class like `Named`, not a second Node base.
- Grep-proof used in CI review: `grep -rn "public  *virtual" core/include/bertini2/function_tree/`
  must stay empty.
