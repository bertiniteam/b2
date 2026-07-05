# The structured output directory: record format `ledgerrec/1`

*The 30-year contract for bertini2's durable records (the structured-output-directory
arc; ADR-0045).  This document plus a directory is a complete, interpretable artifact —
no software required.  A copy of the essentials travels inside every directory as its
`README.txt`.*

## Layout

```
<output directory>/
  README.txt      what this is + this spec's essentials, self-contained
  results.json    the declared results, pretty-printed and self-complete:
                  {"results": {...}, "runs": {...}} -- the final results first,
                  then references (by definition id) to the exact system and
                  configs that produced them
  INDEX.txt       one line per run: when, what, how many paths
  history/        the records: JSONL, one file per writing session, date-named
  definitions/    content-addressed definitions, grouped by kind:
                  definitions/<kind>/<2 hex>/<kind>-<full digest>.<ext>
```

`history/` and `definitions/` are the **source of truth**; the four top-level files are
derived, rebuildable views.  The directory is fully self-contained: deleting it can
never break anything anywhere; the sole consequence is forgetting (worst case:
recomputing).

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
  - **`systems/`** are JSON documents `{"schema", "digest", "encoding", "rendering"}`
    -- one `json.load` away.  The `encoding` is the system's canonical form
    (`b2sysenc/<n>`): exact and complete (every block -- slices, randomization,
    blends -- patch, and gamma survives), and it is the *preimage* of
    `System::ContentDigest()`, so the definition id **equals** the system's identity
    digest and `jq -r .encoding <file> | sha256sum` reproduces it.  The `rendering`
    is classic-style text for eyes; it cannot express the block structure and is
    never an identity (the same text also rides in the run header as
    `target_rendering`).
  - **`configs/`** are JSON with the digest embedded; the digest contract is over the
    `b2cfgenc/<n>` canonical text the JSON is derived from.
  - **`givens/`** are stored byte-exact (fidelity: a CLI input file is *your* file);
    the id is the SHA-256 of the bytes, so plain `sha256sum` reproduces it.

Definitions are written atomically (write-temp, rename) and idempotently (equal content
= equal path, so concurrent writers race benignly; no locks exist or are needed).

## History (`history/`)

Files are named `YYYYMMDD_HHMMSS-pid<pid>[suffix].jsonl` — **one writer per file,
ever** (the name is claimed with exclusive create).  Each line is one JSON object with
a `kind`.  A reader MUST tolerate a torn final line (a crash mid-append) and MUST
treat a torn interior line as corruption.  A reader MUST preserve records whose `kind`
it does not recognize (the operation vocabulary is open; in-progress algorithms mint
new kinds without coordination).

### Record kinds

Every record carries `"schema": "ledgerrec/1"` on its run headers; other records are
scoped by their run.

- **`run`** — one solve/continuation/operation instance:
  `{"kind":"run", "schema":"ledgerrec/1", "when":"YYYY-MM-DD HH:MM", "run":<run id>,
    "op":<operation name, default "solve">, "ask":{...}, "target_object":<definition id>,
    "target_digest":<system content digest -- equals target_object>,
    "target_rendering":<classic-style text, for eyes only, NOT identity>,
    "producer":{"name","version","commit"}, "num_paths":N, ...op-specific fields...}`
  `producer` says which software wrote the record (commit is `unknown` for builds
  outside a git checkout); it is descriptive ONLY -- never part of the ask identity,
  so a newer build answering the same ask recalls rather than recomputes.
  The `ask` is what was requested (target digest + config digest(s) + seed); the run id
  is a hash of the ask.  Op-specific fields describe how start points arise: recorded
  values (`start_points`, exact coordinate text), or a reference to an ancestor run
  (`start_run` + `start_indices`).
- **`track`** — one continued path:
  `{"kind":"track", "run":<run id>, "index":i,
    "status":"success"|"diverged"|"failed",
    "endpoint":[[re,im],...] | null, "start":<provenance>}`
  `status` is the coarse verdict: `diverged` means the path went to infinity — a clean
  `GoingToInfinity` verdict or a deliberate security truncation near infinity
  (`SecurityMaxNormReached`) — an ANSWER, not a failure; `failed` means the tracker
  gave up.  Every tracked path is recorded, whatever its outcome, so a run can be
  audited path-by-path.  The exact codes ride along as
  `pre_endgame_success_code`/`endgame_success_code` (integers) and
  `*_success_code_name` (fixed canonical names — the durable rendering).
  Coordinates are decimal strings at full computed precision (`[real, imaginary]`
  pairs, one per variable, in the target's variable order).  `endpoint` is the
  INTERNAL point (labels: the run header's `variables`); successful/diverged paths
  also carry `endpoint_user`, the dehomogenized point in the USER's coordinates
  (labels: `variables_user`) -- the ones audits should read.  `start` is either
  `{"kind":"start_label","index":i}` (a canonical start-system label — provenance
  bottoms out) or `{"kind":"point_ref","run":<id>,"index":i}` (a chain link into an
  ancestor run's endpoint).
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

### Provenance

Walk a point's ancestry by following `track.start` references backward until a
`start_label` or `given`.  Ids are content-derived, so records from different
directories/sessions merge by file concatenation; matching points ACROSS lineages is
never automatic (float results are not bit-reproducible) — an adopted match is its own
recorded assertion, not a `point_ref`.

## Versioning

`ledgerrec/1` names this format.  Any change to record shapes or file conventions bumps
the version; readers encountering a newer version should read what they understand and
preserve the rest.  The companion identity encodings carry their own versions:
`b2sysenc/<n>` (system canonical encoding, ADR-0042), `b2cfgenc/<n>` (config encoding,
ADR-0043), `b2rand/<n>` (seed-rooted draw derivation, ADR-0044).
