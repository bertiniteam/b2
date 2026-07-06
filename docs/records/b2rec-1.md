# The structured output directory: record format `b2rec/1`

*The 30-year contract for bertini2's durable records (the structured-output-directory
arc; ADR-0045).  This document plus a directory is a complete, interpretable artifact —
no software required.  A copy of the essentials travels inside every directory as its
`README.txt`.*

## Layout

```
<output directory>/
  README.txt      what this is + this spec's essentials, self-contained
  INDEX.txt       one line per run: when, what, how many paths done
  history/        WHAT WAS ASKED, WHEN: run headers, declared results,
                  annotations -- JSONL, one file per writing session, date-named;
                  every line small, referring into the other stores by id
  results/        WHAT WAS COMPUTED: one append-only JSONL file per run,
                  results/<2 hex of run id>/<run id>.jsonl -- a self-description
                  header line, then one line per completed path (endpoint AND
                  its per-path metadata together)
  definitions/    WHAT THINGS ARE: content-addressed inputs, grouped by kind:
                  definitions/<kind>/<2 hex>/<kind>-<full digest>.<ext>
```

`history/`, `results/`, and `definitions/` are the **source of truth** -- three stores,
separated by concern; `README.txt` and `INDEX.txt` are derived, rebuildable views.  The
directory is fully self-contained: deleting it can never break anything anywhere; the
sole consequence is forgetting (worst case: recomputing).

## Definitions (`definitions/`)

One file per definition, filed as
`definitions/<kind>/<first 2 hex of digest>/<kind singular>[-<role>]-<full 64-hex digest>.<ext>`,
e.g. `definitions/systems/03/system-03958a...7f.json` or
`definitions/givens/34/given-cli_input-3468cd...9b.txt`.  The id always sits between
the LAST `-` and the extension; the optional role label (givens carry `cli_input` or
`start_points`) is presentation, never identity.

- The **kind** folder is for the browsing human; bertini writes `systems/`, `configs/`,
  and `givens/`, and the vocabulary is open (any lowercase name).  The kind is
  **presentation, not identity**: records reference bare ids, and an id resolves to the
  unique file whose name carries it, whatever kind it lives under.
- The **filename carries the full digest** -- recovering a definition's id never
  requires string concatenation -- plus the kind and an honest extension (`.json` for
  JSON, `.txt` for text), so a file copied out of the store stays identified.
- The **two-hex shard** folder exists purely so no directory grows unbounded (a
  hundred-thousand-target sweep must not melt `systems/`).
- Every definition is **self-verifying**, per kind:
  - **`systems/`** are JSON documents `{"schema", "digest", "system", "encoding"}`
    -- one `json.load` away.  The `system` value is the structured parts view:
    fields for `variable_groups`, `hom_variable_groups`, `homogenizing_variables`,
    `path_variable`, `named_subexpressions`, `functions` (expression strings), and
    `is_patched`/`num_patches`.  The `encoding` is the system's canonical form
    (`b2sysenc/<n>`): exact and complete (every block -- slices, randomization,
    blends -- patch, and gamma survives), and it is the *preimage* of
    `System::ContentDigest()`, so the definition id **equals** the system's identity
    digest and `jq -r .encoding <file> | sha256sum` reproduces it.  Classic (Bertini
    1) syntax appears nowhere in the records: it is an input/compatibility format
    for replicating results in that software (`bertini.to_classic_input` produces it
    on demand), not an output format -- and it cannot express the block structure.
  - **`configs/`** are JSON with the digest embedded; the digest contract is over the
    `b2cfgenc/<n>` canonical text the JSON is derived from.
  - **`givens/`** are stored byte-exact (fidelity: a CLI input file is *your* file);
    the id is the SHA-256 of the bytes, so plain `sha256sum` reproduces it.

Definitions are written atomically (write-temp, rename) and idempotently (equal content
= equal path, so concurrent writers race benignly; no locks exist or are needed).

## History (`history/`)

The narrative: what was asked, when.  History contains **no per-path lines** -- the
computed paths live in `results/`, referred to by the run header.  Files are named
`YYYYMMDD_HHMMSS-pid<pid>[suffix].jsonl` — **one writer per file, ever** (the name is
claimed with exclusive create).  Each line is one JSON object with a `kind`.  A reader
MUST tolerate a torn final line (a crash mid-append) and MUST treat a torn interior
line as corruption.  A reader MUST preserve records whose `kind` it does not recognize
(the operation vocabulary is open; in-progress algorithms mint new kinds without
coordination).

### Record kinds

Every record carries `"schema": "b2rec/1"` on its run headers; other records are
scoped by their run.

- **`run`** — one solve/continuation/operation instance:
  `{"kind":"run", "schema":"b2rec/1", "when":"YYYY-MM-DD HH:MM", "run":<run id>,
    "op":<operation name, default "solve">, "ask":{...},
    "target_digest":<the solved system's content digest -- also the id of its
                     definitions/systems/ file; dereferencing is the reader's job>,
    "producer":{"name","version","commit"}, "num_paths":N,
    "results_file":"results/<2 hex>/<run id>.jsonl", ...op-specific fields...}`
  `producer` says which software wrote the record (commit is `unknown` for builds
  outside a git checkout); it is descriptive ONLY -- never part of the ask identity,
  so a newer build answering the same ask recalls rather than recomputes.
  The `ask` is what was requested (target digest + homotopy digest + config digest(s)
  + seed); the run id is a hash of the ask.  BOTH the target system and the homotopy
  actually tracked are archived in `definitions/systems/` under their digests -- the
  homotopy's exact coefficients (gamma, start-system constants, blend structure) are
  what the paths followed, and for user-built homotopies the archived encoding is the
  only complete record of them.  The `seed` is the run's own EFFECTIVE seed: a solve
  given no explicit seed derives one at solve start (deterministically chained from
  the session master), so the recorded value reproduces the run standalone -- it is
  never a session master that silently under-determines a mid-session run.  Op-specific fields
  describe how start points arise: recorded values (`start_points`, exact coordinate
  text), or a reference to an ancestor run (`start_run` + `start_indices`).
- **`result`** — a declared deliverable (the signal/noise line):
  `{"kind":"result", "name":<string>, "description":<string>, "when":...,
    "points":[{"run":..,"index":..},...], "value":<any JSON, optional>}`
  Re-declaring a name replaces it (newest wins).  `history/` is everything; only
  declared results are "what I cared about".
- **`annotation`** — metadata attached to a point:
  `{"kind":"annotation", "point":{"run":..,"index":..}, "key":<string>, "value":<JSON>}`
- **`given`** — externally supplied data entering the provenance graph:
  `{"kind":"given", "source":<definition id>, ...}`.  Provenance bottoms out honestly
  at the boundary of what was computed here.

## Results (`results/`)

The payload store: one append-only JSONL file per run at
`results/<2 hex of run id>/<run id>.jsonl`, sharded like `definitions/` so no directory
grows unbounded.  One writer per run in practice (the solve's manager rank); a
concurrent identical ask appends identical lines, which last-wins reading tolerates.
Same torn-line rules as history.  Line 1 is the file's **self-description**:

- **`results_header`** — `{"kind":"results_header", "schema":"b2rec/1",
  "run":<run id>, "ask":{...}}` -- a wandering results file says what run and ask it
  answers, just as a wandering definition self-verifies.

Then one line per completed path, appended as paths finish (a partially-computed run
reads honestly: what is here is done):

- **`path`** — one continued path:
  `{"kind":"path", "run":<run id>, "index":i,
    "status":"success"|"diverged"|"failed",
    "endpoint":[[re,im],...] | null, "start":<provenance>, ...per-path metadata...}`
  `status` is the coarse verdict: `diverged` means the path went to infinity — a clean
  `GoingToInfinity` verdict or a deliberate security truncation near infinity
  (`SecurityMaxNormReached`) — an ANSWER, not a failure; `failed` means the tracker
  gave up.  Every tracked path is recorded, whatever its outcome, so a run can be
  audited path-by-path.  The per-path metadata travels WITH the point (the data and
  its facts are one thing): the exact codes as
  `pre_endgame_success_code`/`endgame_success_code` (integers) and
  `*_success_code_name` (fixed canonical names — the durable rendering), cycle number,
  precision, timings.  Coordinates are decimal strings at full computed precision
  (`[real, imaginary]` pairs, one per variable, in the target's variable order).
  `endpoint` is the INTERNAL point (labels: the run header's `variables`);
  successful/diverged paths also carry `endpoint_user`, the dehomogenized point in the
  USER's coordinates (labels: `variables_user`) -- the ones audits should read.
  `start` is either `{"kind":"start_label","index":i}` (a canonical start-system
  label — provenance bottoms out) or `{"kind":"point_ref","run":<id>,"index":i}` (a
  chain link into an ancestor run's endpoint).

The payload vocabulary is open like the record vocabulary: future operations may
append other payload kinds (e.g. component-membership verdicts) to their run's file.

**Write ordering**: a run's results file (with its header) is created BEFORE the
history run header that refers to it, and each path line is written before anything
refers to that path -- a reference in this format never dangles; at worst an orphaned
payload awaits a reference that never came (harmless, ignorable).

### Provenance

Walk a point's ancestry by following path records' `start` references backward until a
`start_label` or `given`.  Ids are content-derived, so records from different
directories/sessions merge by file concatenation; matching points ACROSS lineages is
never automatic (float results are not bit-reproducible) — an adopted match is its own
recorded assertion, not a `point_ref`.

## Versioning

`b2rec/1` names this format.  Any change to record shapes or file conventions bumps
the version; readers encountering a newer version should read what they understand and
preserve the rest.  The companion identity encodings carry their own versions:
`b2sysenc/<n>` (system canonical encoding, ADR-0042), `b2cfgenc/<n>` (config encoding,
ADR-0043), `b2rand/<n>` (seed-rooted draw derivation, ADR-0044).
