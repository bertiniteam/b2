# ADR-0042: System content identity — canonical encoding, stable digest, Seal(), and interning

**Status:** Accepted
**Date:** 2026-07-02

## Context

Node-level hash-consing shipped in PR #25 (ADR-0027, ADR-0028): every node factory routes through a
process-global intern table keyed by memoized structural hash, so structurally-equal subtrees collapse
to one shared object. ADR-0027's roadmap named the next rung **E4 — intern identical Programs**
("content-hash + weak table; the Program-level analogue of node hash-consing").

The motivation goes beyond memory. The long-term goal is a **database of solutions**: a persistent
provenance ledger in which every continuation records its inputs, operation, and outputs *by id* —
audit is walking ancestors, recomputation is walking descendants, restart is the recorded status.
That architecture stands on identity: provenance must *reference* shared, interned definitions by
handle, never embed copies — the interning **is** the identity mechanism. A lineage id for an endpoint
is its start-label transformed by the chain of homotopy ids, which requires each homotopy — including
its recorded randomness (gamma, patch and slice coefficients) — to have a canonical, stable identity.

Two facts about the existing machinery made this non-trivial:

1. **The node hash is in-process only.** `Node::Hash()` mixes `typeid(...).hash_code()` and
   `std::hash<std::string>` — stable within one process run, *not* across runs, standard libraries, or
   compilers. Perfect for the intern table; unusable as a persistent key.
2. **Systems are mutable.** Nodes could hash-cons at construction because they are immutable
   (ADR-0011). A `System` accretes state over many calls (`AddFunction`, `Homogenize`, `AutoPatch`,
   variable grouping…) — constructing a system by mutation is the natural authoring style and stays.
   Identity therefore cannot attach at construction.

## Decision

One canonical substrate, two identity layers, plus seal semantics and interning at both the System and
Program levels.

### The substrate: a purpose-built, versioned canonical encoding

`CanonicalEncoding(System)` produces an exact, deterministic text encoding; `System::ContentDigest()`
is its SHA-256. This digest is the persistent identity key (the lineage-id building block).

Rejected substrates, and why:

- **Boost.Serialization archive bytes** — archives embed the boost archive-format version and
  per-class versions (bytes rotate on a Boost upgrade), and `shared_ptr` tracking ids depend on
  serialization *order and history*, so equal content can serialize to different bytes.
- **The classic printer (ADR-0012)** — printed forms are documented as a cosmetic, test-pinned
  surface; a printer tweak would silently rotate every persistent digest. The printer guarantees
  reparse-*value*-invariance, not a bit-exact spec.
- **`boost::uuids::detail::sha1`** — the `detail::` namespace carries no stability contract across
  Boost versions, and SHA-1 is the weaker choice for a forever-key.

The encoding is S-expression-like with fixed string kind tags (never `typeid`), exact coefficient
strings (`mpz`/`mpq` `.str()`; mpfr full-precision digits plus the stored precision), and
first-occurrence back-references assigned in traversal order — so equal trees encode identically
regardless of incidental pointer sharing, while heavily shared DAGs stay linear-size. The System
encoding starts with a **format-version header** (`b2sysenc/1`) and the identity-affecting
session-global canonicalization settings (monomial order, canonicalize-by-default, power-folding):
digests produced under different settings correctly compare unequal rather than silently colliding.
Any change to the encoding spec **must** bump the version and update the golden-digest fixture in the
same commit — digest drift is the one unrecoverable failure mode for a database keyed on digests.

**Everything evaluation-relevant is identity**: functions, variable groups with their types and
time-order, path variable, parameters, patch coefficients, randomization matrices, blend coefficients
(gamma). Randomness is *recorded as part of identity*, not excluded — two homotopies with different
gammas are different homotopies. Pre-homogenization functions are included (observable through
user-coordinate Jacobians). Excluded is only transient/derived state: current variable values,
working precision, differentiation caches, SLPs.

Variable **names are load-bearing** (variables intern globally by name); identity is not up to
alpha-equivalence.

### In-memory identity is derived from the substrate

`System::Hash()` folds the first 8 digest bytes; `System::IsSame()` compares digests. The invariant
`IsSame ⇒ equal Hash` holds trivially, and there is no second hand-written structural comparator over
five block types to keep in sync. Systems are compared at seal/intern time, never in a hot loop, and
sealing memoizes the digest.

### Seal semantics: hashcons-on-freeze

`System::Seal()` memoizes the digest and flips the system into a sealed state; structural mutators
throw `std::logic_error` thereafter. `ContentDigest()` works on any system (computed fresh when
unsealed). Transient operations — `SetVariables`, `precision(unsigned)`, `Differentiate()`,
evaluation — remain allowed on a sealed system. Copying a sealed system yields an **unsealed** copy
(the copy-on-write escape hatch).

The name is "Seal", not "Freeze": ADR-0027's E3 "freeze set" is evaluation-time input currying, a
different concept. They compose — sealing forbids only structural mutation.

A `FrozenSystem` wrapper type was rejected: it would fork the entire System-consuming API surface
(trackers, endgames, blocks holding operand systems) for no semantic gain.

### Interning

- **Systems:** `InternSystem(candidate)` returns the live sealed representative with equal
  `ContentDigest()` if one exists, else seals and registers the candidate. The table is a
  process-global mutex-guarded `map<Digest256, weak_ptr<const System>>`, expired entries pruned on
  touch — mirroring `node::Intern`. The full 256-bit key needs no disambiguation chain.
- **Programs (E4 proper):** `SLPProgram` gains an in-process `ContentHash()`/`SameContent()` and a
  weak intern table wired into `SLPCompiler::Compile` and deserialization, so identical compiled tapes
  are shared. Program identity is compiler-version-dependent, hence in-memory only — the System
  digest is the persistent key. `SLPMemory` stays per-facade (the ADR-0027 threading contract).

### Re-interning on load

Deserialization deliberately bypasses `Intern` (nodes are constructed raw). Left alone, that defeats
hash-consing across sessions: on a System-table miss the loaded system would become a representative
whose *nodes* live outside the node intern table, so nothing built afterward — and nothing loaded from
a second archive — would share with it. Therefore:

- **Node level:** `node::Reintern(root)` rebuilds the loaded DAG post-order through the interning
  factories, with a per-class rebuild hook and a pointer-keyed memo that preserves the archive's
  internal sharing (linear pass). Rebuilt nodes acquire fresh memoized state. This pass runs only at
  deserialization time; the construction-internal `Intern` bypass is untouched.
- **System level:** after load, every node-holding member (block functions, blend coefficients,
  operand systems recursively, variable groups, homogenizing variables, path variable, explicit
  parameters) is remapped through the same memo, then the system goes through `InternSystem`.
  `LoadSystemUnified` packages load → reintern → intern.

## Consequences

- A System (and hence a homotopy, randomness included) has a stable, cross-session, content-derived
  identity — the L2 lineage-id building block for the provenance ledger — plus cheap in-memory
  `Hash`/`IsSame`.
- Identical compiled Programs are shared process-wide (memory win under parallel tracking), and
  loaded artifacts unify with live ones instead of forking the intern universe.
- The digest is only as stable as the encoding spec; the golden-digest fixture test
  (`core/test/classes/data/system_digest_fixture.txt`) fails loudly on any drift, and the versioned
  header turns an intentional spec change into a new keyspace rather than silent corruption.  The
  bump itself is also under test: the append-only version registry
  (`core/test/classes/data/system_encoding_versions.txt`) maps every `b2sysenc` version ever used to
  a hash of the fixture recipes' encoding texts (version token excluded), so regenerating the
  fixture without bumping the version fails the registry test.
- Sealing is opt-in: existing mutable workflows are untouched, and factory outputs
  (`MakeHomotopy`, `Randomize`, zero-dim internals) are not auto-sealed.
- Value-equal coefficients stored at different precisions are *distinct* identities (their
  high-precision evaluations differ).

## Still open (follow-ups, not done here)

- Auto-sealing factory outputs once the ecosystem is comfortable with seal semantics.
- Alpha-equivalent identity (positional variable encoding) if renamed-but-equal systems ever need to
  unify — a spec version bump.
- The provenance ledger itself: this ADR delivers the identity mechanism it stands on, not the ledger.
