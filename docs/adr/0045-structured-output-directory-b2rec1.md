# ADR-0045: The structured output directory — record schema b2rec/1 and the C++ OutputDirectory

**Status:** Accepted
**Date:** 2026-07-03

> **Update (2026-07-06, rung-7 shakedown):** the layout evolved from *two* truth stores
> (`history/` + `definitions/`, with a single derived `results.json`) to **three**:
> `history/` (what was asked, when — run headers, declared results, annotations),
> `results/` (one append-only JSONL file per run — every computed path with its endpoint
> and metadata), and `definitions/` (content-addressed inputs).  `results.json` and
> `RESULTS.txt` are **retired** — `history/` no longer carries per-path lines, so the
> "database" concern moved out of the journal into `results/`.  `README.txt` + `INDEX.txt`
> remain the only derived views.  The current spec is `docs/records/b2rec-1.md`; the
> rationale is in the arc doc (`arcs/structured-output-directory.md`, rounds 7–10).  The
> rest of this ADR describes the original two-store design.

## Context

The structured-output-directory arc (rung 3) needs its storage layer in core C++, so
both faces of bertini2 — the CLI and Python — write the same durable records.  The
requirements were set by the arc's design sessions: interactable WITHOUT special
software (plain text is the source of truth; anything fancier is a derived index),
self-documenting for a decade (a zipped directory on Zenodo must explain itself),
humane above the fold ("a mathematician is smart but lazy"), safe under crashes and
concurrent writers on parallel filesystems, and an OPEN operation vocabulary so
in-progress algorithms (NID regeneration) can mint record kinds without coordination.
The pure-Python pilot (`prototypes/ledger_v0`) validated every choice; this ADR makes
them a specification and a core implementation.

## Decision

- **`b2rec/1`** is the record format, specified in `docs/records/b2rec-1.md`
  and carried inside every directory as its `README.txt`.  Layout: `history/`
  (append-only JSONL, one date-named file per writing session, ONE writer per file,
  claimed by exclusive create) and `definitions/` (content-addressed under
  human-navigable kind folders, filenames carrying kind + full digest + honest
  extension — `systems/` are JSON embedding the canonical encoding, the digest
  preimage, so the definition id IS the system's identity digest; atomic
  write-temp+rename, idempotent — concurrent writers race benignly, no locks exist)
  are the source of truth; `README.txt`, `INDEX.txt`, `results.json`
  are derived, rebuildable views.  Record kinds: `run`, `track`, `result`,
  `annotation`, `given`; readers MUST tolerate a torn final line, MUST treat torn
  interior lines as corruption, and MUST preserve unknown kinds.
- **`bertini::records::OutputDirectory`**
  (`core/include/bertini2/records/output_directory.hpp`) implements it: definitions
  put/get, session-journal append (flush per record: a walltime kill loses at most the
  line in flight), tolerant scan, and the derived-view renderers (INDEX; results.json
  for code AND eyes: pretty-printed, with a "runs" section referring to the
  system/config definitions each result came from.  RESULTS.txt was retired once
  results.json became human-readable -- one results file, no duplication).
- **Boost.JSON, header-only**: `<boost/json/src.hpp>` is included in exactly one TU,
  so no new link component is required on any platform (CI builds Boost from source
  with a fixed component list).
- **Cross-implementation compatibility is tested, both directions**: the C++ suite
  writes an example directory the Python pilot reads (endpoint decode, sha256
  self-verification, results.json load), and reads one the pilot writes.  The byte
  formats agree where it matters — definitions are byte-identical by construction
  (content-addressing is the same hash of the same bytes); records are
  JSON-value-equivalent (serializer whitespace/float-rendering differences between
  boost::json and Python's json are immaterial to a format whose numbers-of-record are
  strings).

## Consequences

- Rung 4 (the solver emission seam) has its storage layer: solvers append records and
  consult the directory with no further storage design.
- The CLI inherits the same directory the Python face writes — one format, two faces.
- The derived views are quality-of-life, not truth: any of the four top-level files can
  be deleted and regenerated from `history/` + `definitions/`.
- The cross-impl bridge tests run fully only where both implementations are present
  (local dev, C++ CI stages); the read-the-python-example leg self-skips in wheel CI.
- The v0 pilot is now redundant at the storage layer; it remains the executable spec
  for the memoized-solve semantics until rung 5 retires it entirely.
