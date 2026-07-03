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
                  then references (by definition id) to the system rendering and
                  configs that produced them
  INDEX.txt       one line per run: when, what, how many paths
  history/        the records: JSONL, one file per writing session, date-named
  definitions/    content-addressed definitions: definitions/<2 hex>/<rest of digest>
```

`history/` and `definitions/` are the **source of truth**; the four top-level files are
derived, rebuildable views.  The directory is fully self-contained: deleting it can
never break anything anywhere; the sole consequence is forgetting (worst case:
recomputing).

## Definitions (`definitions/`)

One file per definition, named by SHA-256.  Two id conventions:

- **self-verifying**: the id is the SHA-256 of the file's bytes (`sha256sum` checks it);
- **external**: the id is the object's own content digest (e.g. a System's
  `ContentDigest()`, whose preimage is its canonical encoding `b2sysenc/<n>`), and the
  file holds a human-readable rendering (classic input).

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
    "num_paths":N, ...op-specific fields...}`
  The `ask` is what was requested (target digest + config digest(s) + seed); the run id
  is a hash of the ask.  Op-specific fields describe how start points arise: recorded
  values (`start_points`, exact coordinate text), or a reference to an ancestor run
  (`start_run` + `start_indices`).
- **`track`** — one continued path:
  `{"kind":"track", "run":<run id>, "index":i, "status":"success"|"failed",
    "endpoint":[[re,im],...] | null, "start":<provenance>}`
  Coordinates are decimal strings at full computed precision (`[real, imaginary]`
  pairs, one per variable, in the target's variable order).  `start` is either
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
