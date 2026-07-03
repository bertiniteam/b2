# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""Ledger v0: a content-addressed definition store + append-only JSONL history.

Plain files are the source of truth (arc: structured-output-directory).  Layout is for
humans first ("a mathematician is smart but lazy"): `history/` holds date-named record
files, `definitions/` holds content-addressed definitions, README.txt explains the
directory in place, and INDEX.txt summarizes every run -- `cat INDEX.txt` answers
"what's in here?".

Definitions are written atomically (tmp + rename) and idempotently (same content = same
path, so races are benign).  History files are one-writer-per-file; the reader tolerates
a torn final line (the crash-mid-append case) and refuses torn lines anywhere else.
"""

import hashlib
import json
import os
import time
import uuid
from pathlib import Path

SCHEMA = "ledgerrec/0"


def _sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class Ledger:
    """One structured output directory: definitions/ + history/ + README.txt + INDEX.txt."""

    README = """This directory is a structured output directory: the durable, self-contained
record of numerical algebraic geometry computations (polynomial-system solves by
homotopy continuation).  Written by bertini2 (ledger v0 prototype, record schema
ledgerrec/0).  It needs no software to read, and you are free to delete it -- the
only consequence is recomputing.

LAYOUT
  RESULTS.txt   the declared results, for EYES (never parse this).
  results.json  the same results, for CODE: one json.load away.
                Most readers start AND END with these two.
  INDEX.txt     one line per run: when, what was solved, how many paths.
  history/      the records: JSON, one object per line (JSONL), one file per writing
                session, named by date.  Read with eyes, grep, jq, or
                pandas.read_json(..., lines=True).
  definitions/  the things records refer to (polynomial systems, etc.), stored once
                each in a file named by the SHA-256 of its content -- so records can
                reference them exactly, and `sha256sum` verifies them.  Mostly
                human-readable input files.

HOW REFERENCES WORK (one direction, no cycles)
  INDEX.txt names runs  ->  runs and their per-path `track` lines live in history/
  ->  track lines reference ancestor points by {run, index} (that is the provenance
  graph)  ->  and every record references its systems by hash into definitions/.
  RESULTS.txt / results.json are DERIVED views of the above, rebuildable at will.

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
"""

    def __init__(self, root):
        self.root = Path(root)
        (self.root / "definitions").mkdir(parents=True, exist_ok=True)
        (self.root / "history").mkdir(parents=True, exist_ok=True)
        readme = self.root / "README.txt"
        if not readme.exists():
            readme.write_text(self.README)

    # ---- object store -------------------------------------------------------------

    def _object_path(self, object_id: str) -> Path:
        return self.root / "definitions" / object_id[:2] / object_id[2:]

    def put_object(self, data, object_id=None) -> str:
        """Store bytes/str content-addressed; returns the object id.

        If object_id is given (e.g. a System's content_digest, whose preimage is the
        canonical encoding rather than these bytes), the file is stored under that id;
        otherwise the id is the sha256 of the bytes (self-verifying).
        """
        if isinstance(data, str):
            data = data.encode("utf-8")
        oid = object_id if object_id is not None else _sha256_hex(data)
        path = self._object_path(oid)
        if path.exists():
            return oid  # idempotent: content-addressed writes never conflict
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp-%s" % uuid.uuid4().hex[:8])
        tmp.write_bytes(data)
        os.replace(tmp, path)  # atomic on POSIX
        return oid

    def get_object(self, object_id: str) -> bytes:
        return self._object_path(object_id).read_bytes()

    def has_object(self, object_id: str) -> bool:
        return self._object_path(object_id).exists()

    # ---- journals -----------------------------------------------------------------

    def append(self, record: dict):
        """Append one record to THIS session's journal (opened lazily, one per Ledger
        instance / process session).  One writer per file, ever; the date-stamped name
        makes `ls history/` read as a history, not confetti."""
        if getattr(self, "_journal", None) is None:
            stamp = time.strftime("%Y%m%d_%H%M%S")
            base = self.root / "history"
            for suffix in [""] + ["%c" % c for c in range(ord("b"), ord("z"))]:
                path = base / ("%s-pid%d%s.jsonl" % (stamp, os.getpid(), suffix))
                try:
                    path.touch(exist_ok=False)   # claim the name exclusively
                    break
                except FileExistsError:
                    continue
            self._journal = Journal(path)
        self._journal.append(record)
        if record.get("kind") == "run":
            self.refresh_index()

    def describe(self) -> str:
        """One human line: how much is here.  Printed by demos at exit so the records'
        location is never a mystery."""
        self.refresh_index()
        history = list((self.root / "history").glob("*.jsonl"))
        n_defs = sum(1 for _ in (self.root / "definitions").glob("*/*"))
        n_records = len(self.scan())
        return "%d records in %d history file(s), %d definition(s), at %s" % (
            n_records, len(history), n_defs, self.root.resolve())

    def refresh_index(self):
        """(Re)write INDEX.txt: one line per run -- when, op, what, path counts.  Derived
        from the records (rebuildable at will); the lazy mathematician's front door."""
        records = self.scan()
        runs = [r for r in records if r.get("kind") == "run"]
        track_counts = {}
        for r in records:
            if r.get("kind") == "track":
                counts = track_counts.setdefault(r["run"], {"success": 0, "failed": 0})
                counts[r["status"]] = counts.get(r["status"], 0) + 1
        lines = ["what has been solved here (newest last; details in history/):", ""]
        for run in runs:
            counts = track_counts.get(run.get("run"), {})
            n_done = sum(counts.values())
            status = "%d/%s paths done" % (n_done, run.get("num_paths", "?"))
            if counts.get("failed"):
                status += " (%d failed)" % counts["failed"]
            lines.append("%s  %s  %-8s  %s   [run %s]" % (
                run.get("when", "????-??-?? ??:??"),
                status.rjust(22),
                run.get("op", "solve"),
                self._describe_target(run.get("target_object")),
                run.get("run", "?")))
        (self.root / "INDEX.txt").write_text("\n".join(lines) + "\n")

    def refresh_results(self):
        """(Re)write the two derived result views -- 'what were my solutions?':

        - results.json  for CODE: one json.load away (coordinates keyed by variable
          name, annotations, provenance refs, inline saved values).
        - RESULTS.txt   for EYES only: pretty, never to be parsed.

        Both derived and rebuildable, like INDEX.txt; together the humane analogue of
        bertini1's main_data."""
        records = self.scan()
        runs = {r["run"]: r for r in records if r.get("kind") == "run" and "run" in r}
        tracks = {(r["run"], r["index"]): r for r in records if r.get("kind") == "track"}
        notes = {}
        for r in records:
            if r.get("kind") == "annotation":
                key = (r["point"]["run"], r["point"]["index"])
                notes.setdefault(key, {})[r["key"]] = r["value"]
        declared = {}
        for r in records:
            if r.get("kind") == "result":
                declared[r["name"]] = r        # newest declaration of a name wins

        machine = {}
        lines = ["the results declared here (for eyes only -- code reads results.json;",
                 "full precision + provenance in history/):"]
        for name, rec in declared.items():
            entry = {"declared": rec.get("when"), "description": rec.get("description", "")}
            lines += ["", "== %s   (declared %s) ==" % (name, rec.get("when", "?"))]
            if rec.get("description"):
                lines.append("   %s" % rec["description"])

            if "value" in rec:                 # an inline save(): arbitrary JSON-able thing
                entry["value"] = rec["value"]
                lines.append("  value: %s" % json.dumps(rec["value"]))

            points = []
            for n, ref in enumerate(rec.get("points", []), 1):
                key = (ref["run"], ref["index"])
                track = tracks.get(key)
                if track is None or track.get("status") != "success":
                    lines.append("  point %d: (not computed / failed)" % n)
                    points.append({"provenance": ref, "status": "missing"})
                    continue
                note = notes.get(key, {})
                var_names = self._variable_names(runs.get(ref["run"], {}).get("target_object"))
                coords = {}
                for k, (re_str, im_str, *_prec) in enumerate(track["endpoint"]):
                    var = var_names[k] if k < len(var_names) else "x%d" % k
                    coords[var] = [re_str, im_str]
                points.append({"coordinates": coords, "annotations": note,
                               "provenance": ref, "status": "success"})

                suffix = ("   " + ", ".join("%s=%s" % kv for kv in sorted(note.items()))
                          if note else "")
                lines.append("  point %d  [run %s #%d]%s" % (n, ref["run"], ref["index"], suffix))
                for var, (re_str, im_str) in coords.items():
                    lines.append("    %s = %.12s + %.12s i" % (var, re_str, im_str))
            if points:
                entry["points"] = points
            machine[name] = entry
        if not declared:
            lines.append("(none declared yet)")
        (self.root / "results.json").write_text(json.dumps(machine, indent=1) + "\n")
        (self.root / "RESULTS.txt").write_text("\n".join(lines) + "\n")

    def _variable_names(self, object_id) -> list:
        """Variable names from a stored classic input, for labeling coordinates."""
        if not object_id or not self.has_object(object_id):
            return []
        text = self.get_object(object_id).decode("utf-8", "replace")
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("variable_group") or stripped.startswith("variable "):
                names = stripped.split(None, 1)[1].rstrip(";")
                return [v.strip() for v in names.split(",")]
        return []

    def _describe_target(self, object_id):
        """A one-glance description of a run's target: its function lines, abbreviated."""
        if not object_id or not self.has_object(object_id):
            return "(unknown target)"
        text = self.get_object(object_id).decode("utf-8", "replace")
        functions = [ln.strip().rstrip(";") for ln in text.splitlines() if "=" in ln]
        summary = "; ".join(functions)
        return (summary[:57] + "...") if len(summary) > 60 else summary

    def scan(self):
        """Read every record from every journal.  Torn final lines are skipped (the
        crash-mid-append case); a torn line anywhere else raises (real corruption)."""
        records = []
        for path in sorted((self.root / "history").glob("*.jsonl")):
            lines = path.read_text(encoding="utf-8").splitlines()
            for lineno, line in enumerate(lines):
                if not line.strip():
                    continue
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    if lineno == len(lines) - 1:
                        continue  # torn tail: the write died mid-line; replay ignores it
                    raise ValueError("corrupt journal %s at line %d" % (path, lineno + 1))
        return records

    # ---- record-level queries -----------------------------------------------------

    def find_run(self, ask: dict):
        """The most recent run header matching this ask, or None."""
        found = None
        for rec in self.scan():
            if rec.get("kind") == "run" and rec.get("ask") == ask:
                found = rec
        return found

    def completed_paths(self, run_id: str) -> dict:
        """index -> track record, for every recorded path of the run."""
        done = {}
        for rec in self.scan():
            if rec.get("kind") == "track" and rec.get("run") == run_id:
                done[rec["index"]] = rec
        return done


class Journal:
    """An append-only JSONL writer.  One line per record; flushed per append so a kill
    loses at most the line being written (which scan() tolerates)."""

    def __init__(self, path: Path):
        self.path = path
        self._file = open(path, "a", encoding="utf-8")

    def append(self, record: dict):
        self._file.write(json.dumps(record, separators=(",", ":")) + "\n")
        self._file.flush()

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
