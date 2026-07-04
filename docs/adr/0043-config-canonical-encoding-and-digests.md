# ADR-0043: Configuration structs get canonical encodings and persistent digests

**Status:** Accepted
**Date:** 2026-07-03

## Context

The structured-output-directory arc (`arcs/structured-output-directory.md`; rung 1)
records every computation with its complete defining identity, so that equal asks are
memoizable (resume, cross-run reuse) and every recorded point's provenance is exactly
replayable.  A computation's ask identity is `(operation, target-system digest,
CONFIG digest, seed)` — the system half shipped in ADR-0042; this ADR delivers the
config half.

The same stability problem ADR-0042 solved recurs: a digest that must mean the same
thing across runs, compilers, machines, and years cannot be built from `std::hash`,
`typeid`, numeric enum values, or formatted floating point.  Configs add one new trap of
their own: **doubles**.  Any decimal rendering is a rounding policy; two configs
differing by one ulp must digest differently, and the same config must digest
identically under every C++ runtime's formatting quirks.

## Decision

A new `bertini::records` namespace (the future home of the output-directory machinery)
with one canonical encoder per configuration struct and digests over them
(`core/include/bertini2/records/config_encoding.hpp`):

- **Versioned header**: every digest is SHA-256 over `b2cfgenc/1` + the encoding.  Any
  change to what an encoder emits bumps the version and the golden fixture in the same
  commit — an intentional new keyspace, never silent drift.
- **Exact scalars**:
  - `double` → `d64:<16 hex>`: the IEEE-754 bit pattern, the only single-valued
    rendering of a double (decimal and `%a` hexfloat forms are implementation-variable);
  - `mpq_rational` → exact `.str()`;
  - enums → fixed string tables (`Predictor::RKF45` → `"RKF45"`), never numeric values;
  - strings → netstrings (matching the node encoder's convention).
- **One encoder per struct, fields by name in declaration order**, for all
  tracking-relevant configs: Stepping, Newton, FixedPrecision,
  AdaptiveMultiplePrecision, Predictor; Security, Endgame, PowerSeries, Cauchy,
  TrackBack; Tolerances, MidPath, AutoRetrack, Sharpening, Regeneration,
  PostProcessing, ZeroDim, Meta, EndgameChoice.
- `ConfigDigest(config)` digests one struct; `SettingsDigest(c1, c2, ...)` digests a
  solver's full settings as the ordered concatenation — **order is part of the
  contract** (each solver kind composes in one fixed, documented order).

Deliberate exclusions:

- **`RandomConfig` is not encodable**: the seed is its own slot in the ask identity,
  beside the config digest, never inside it (seed-rooted randomness, arc rung 2).
- **`ZeroDimConfig::num_threads` is excluded**: thread count is transient — it must not
  change what was computed, so it must not change the identity (tested).

## Consequences

- Any solver's full settings reduce to one stable digest — the config half of ask
  identity — pinned by the golden fixture
  `core/test/classes/data/config_digest_fixture.txt` (eight recipes; drift fails
  loudly with bump-the-version instructions).
- **Adding a field to a config struct now requires extending its encoder.**  The
  encoding is a hand-maintained mirror of the struct (like `serialize`); the fixture
  catches encoder-vs-struct drift only when defaults change, so review discipline
  matters.  A reflection-based generator was rejected: C++17 has no reflection, and a
  macro DSL would obscure the one property that matters — that a human can read exactly
  what is identity.
- Pre-existing wart, surfaced not fixed: `SharpeningConfig::sharpendigits` and the
  three `RegenerationConfig` slice tolerances are default-UNINITIALIZED, and
  `SharpeningConfig::function_residual_tolerance` defaults precision-dependently.
  Encoders read whatever is there; fixtures set every field explicitly.  Initializing
  those defaults is separate cleanup.
- MPI verified (per-rung gate): CLI mpirun smoke (counts match serial, both
  start-system families) and the mpi4py zero-dim suite both green with this change.

## Addendum (2026-07-03): the archived form is JSON; the digest preimage is unchanged

Output directories archive a solve's configs under the settings digest as **pretty
JSON** (`records::ConfigTextAsJson` — `{"schema", "digest", "configs": {Name:
{field: value}}}`), derived mechanically from the canonical text with exact values
preserved as strings.  The DIGEST contract stays on the `b2cfgenc/<n>` canonical text
exactly as decided above; the JSON is a view, chosen because the directory's human
and tool surfaces already speak JSON.

Also: the settings digest for a zero-dim/homotopy solve now composes ALL configs the
solve reads (algorithm + tracker + endgame; see ADR-0046) — the per-config encoders
and the golden fixture are untouched by that composition change.

Coordination note: PR #70 adds `CauchyConfig::num_pole_growth_rounds_before_truncation`;
when it merges alongside the records arc, the Cauchy encoder must gain the field and
the fixture regenerates (adding a config field REQUIRES extending its encoder).
