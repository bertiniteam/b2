This directory is a structured output directory: the durable, self-contained
record of numerical algebraic geometry computations (polynomial-system solves by
homotopy continuation).  Written by bertini2 (ledger v0 prototype, record schema
ledgerrec/0).  It needs no software to read, and you are free to delete it -- the
only consequence is recomputing.

LAYOUT
  RESULTS.txt   the declared results, with coordinates.  Most readers start AND END here.
  INDEX.txt     one line per run: when, what was solved, how many paths.
  history/      the records: JSON, one object per line (JSONL), one file per writing
                session, named by date.  Read with eyes, grep, jq, or
                pandas.read_json(..., lines=True).
  definitions/  the things records refer to (polynomial systems, etc.), stored once
                each in a file named by the SHA-256 of its content -- so records can
                reference them exactly, and `sha256sum` verifies them.  Mostly
                human-readable input files.

RECORD FORMAT (schema ledgerrec/0) -- every history line is one JSON object:
  kind="run"        a solve: `ask` (what was requested: target system digest + config),
                    `run` (this run's id), `when`, `num_paths`, and how start points
                    arise (recorded values, or a reference to an ancestor run).
  kind="track"      one continued path: `run`, `index`, `status`, `endpoint`
                    (coordinates as [real, imaginary] decimal-string pairs, full
                    precision), and `start` (its provenance: a start_label, or a
                    point_ref {run, index} into an ancestor run's endpoint).
  kind="annotation" metadata attached to a point: `point` {run, index}, `key`, `value`.
  kind="result"     the declared DELIVERABLES: `name`, `points` [{run, index}, ...].
                    Everything else in history/ is scaffolding; RESULTS.txt renders
                    these with coordinates -- "what were my solutions?".
Chains of runs are walkable: follow track records' `start` references backward until
a start_label -- that is the complete provenance of any point recorded here.
