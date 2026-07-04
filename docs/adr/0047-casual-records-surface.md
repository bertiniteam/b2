# ADR-0047: The casual records surface — solve / save / load, and records as ordinary output

**Status:** Accepted
**Date:** 2026-07-03

## Context

Rung 5 of the structured-output-directory arc: put humane faces on the ADR-0046 seam.
The governing UX decisions were made in the arc sessions: the casual API is verbs with
zero nouns; records are simply the program's output (like a compiler's .o or LaTeX's
.log) — produced always, freely deletable, NO opt-in/opt-out flag; a solution is a
point that carries its metadata; naming must be self-evident (`saved()` was rejected as
unclear — the pair is **save/load**).

## Decision

- **Python** (`bertini.records`, re-exported at top level):
  - `solve(system, seed=None, directory=None, ...)` — builds the default ZeroDim
    solver, attaches the ambient directory, solves (ensure-answered: hydration via the
    seam), auto-declares its finite solutions as a result, and returns a `SolveResult`:
    a claim ticket (run id + directory + solutions) that is safe to drop.
  - `save(thing)` / `save(name, thing)` — a solve result is declared with full
    provenance ({run, index} refs); any JSON-able value is recorded inline; nameless
    saves auto-name by timestamp.  **`load(name)` / `load()`** reads results.json —
    works in any later session, and the same file needs no bertini at all.
  - `Solution(np.ndarray)` — coordinates that remember: `.provenance`/`.annotations`
    ride invisibly; arithmetic yields plain derived points (provenance honestly
    absent).  `records_dir()` gets/sets the ambient directory
    (`BERTINI_RECORDS_DIR`, else `./bertini_output`).
  - Bindings: solver-level `record_to(path)` / `records_run_id()` /
    `num_paths_hydrated()` / `records_path()` / `refresh_results()` on every solver
    (ZDVisitor), and a minimal `_pybertini.records.OutputDirectory` for WRITING through
    the single C++ implementation (reading is plain json — the point of the format).
- **CLI**: records are on by default — `bertini2` writes `bertini_output` beside the
  b1-compatible files, `BERTINI_RECORDS_DIR` overrides, manager rank only, via a new
  `AnyZeroDim::RecordToPath` virtual (default no-op).  A killed run finishes when the
  same command is run again; zero new flags.
- **Tutorial**: `python/docs/source/tutorials/your_records/` (doctest-wired).
- Existing power-user solver objects record only when asked (explicit `record_to` /
  env var): the 600-test suite stays clean, and "always" applies to the two casual
  faces where the no-flag principle lives.

## Consequences

- The casual contract is complete: `sols = pb.solve(sys, seed=42); pb.save(sols)` —
  rerunning is instant, crashes resume, seeds are shareable identities.
- Endpoints render in INTERNAL coordinates (homogenized) in results.json — correct for
  restart, but user-coordinate rendering of results is a noted follow-up (rung 6
  polish), as is `annotate()` and chained `point_ref` provenance (rung 6 proper).
- The Python prototype (`prototypes/ledger_v0`) is superseded for solve/save/load and
  resume; it is retained until rung 6 retires its chain/curve demos.

## Addendum (2026-07-03/04): the surface as it shipped

Grown since the original decision, same principles:

- **Chains**: `solve(B, homotopy=H, start=r1)` — a prior result's solutions chain
  with `point_ref` provenance; raw arrays are archived as a *given* (`given_ref`
  starts); the start-data identity joins the ask.
- **Verbs**: `annotate(point, key, value)`; `solutions_of(run)` (cold-read any run's
  endpoints — CLI-written included — as chainable Solutions); `provenance(point)`
  (the walk back to a start label or given).
- **Navigation/viz**: `runs()`/`tracks()` (pandas DataFrames; coordinates excluded
  from tracks unless asked — million-path scale guard), `provenance_graph()`
  (networkx DiGraph), `plot_chain()` (left-to-right lineage view; aggregates to
  run-level above `max_paths_drawn`).
- **The off switch is one line**: `bertini.recording(False)`; for the CLI/ambient,
  an EMPTY `BERTINI_RECORDS_DIR`.  Off = bare solve, provenance honestly absent.
- **Track verdicts are three-way** (success / diverged / failed — truncation is a
  verdict, not a failure) with SuccessCode names recorded beside the integers;
  tracks carry `endpoint_user` (+ run headers `variables_user`) so audits read the
  user's coordinates.  `results.json` is the single results file: pretty-printed,
  self-complete (`{"results", "runs"}` with definition refs); RESULTS.txt retired.
- Run headers carry `producer {name, version, commit}` — descriptive, never identity.
