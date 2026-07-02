# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""Ledger v0: a content-addressed object store + append-only JSONL journals.

Plain files are the source of truth (arc: structured-output-directory).  Objects are
written atomically (tmp + rename) and idempotently (same content = same path, so races
are benign).  Journals are one-writer-per-file; the reader tolerates a torn final line
(the crash-mid-append case) and refuses torn lines anywhere else.
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
    """One ledger directory: objects/ (definitions) + journals/ (records)."""

    def __init__(self, root):
        self.root = Path(root)
        (self.root / "objects").mkdir(parents=True, exist_ok=True)
        (self.root / "journals").mkdir(parents=True, exist_ok=True)

    # ---- object store -------------------------------------------------------------

    def _object_path(self, object_id: str) -> Path:
        return self.root / "objects" / object_id[:2] / object_id[2:]

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
        makes `ls journals/` read as a history, not confetti."""
        if getattr(self, "_journal", None) is None:
            stamp = time.strftime("%Y%m%d_%H%M%S")
            base = self.root / "journals"
            for suffix in [""] + ["%c" % c for c in range(ord("b"), ord("z"))]:
                path = base / ("%s-pid%d%s.jsonl" % (stamp, os.getpid(), suffix))
                try:
                    path.touch(exist_ok=False)   # claim the name exclusively
                    break
                except FileExistsError:
                    continue
            self._journal = Journal(path)
        self._journal.append(record)

    def describe(self) -> str:
        """One human line: how much is here.  Printed by demos at exit so the records'
        location is never a mystery."""
        journals = list((self.root / "journals").glob("*.jsonl"))
        n_objects = sum(1 for _ in (self.root / "objects").glob("*/*"))
        n_records = len(self.scan())
        return "%d records in %d journal file(s), %d object(s), at %s" % (
            n_records, len(journals), n_objects, self.root.resolve())

    def scan(self):
        """Read every record from every journal.  Torn final lines are skipped (the
        crash-mid-append case); a torn line anywhere else raises (real corruption)."""
        records = []
        for path in sorted((self.root / "journals").glob("*.jsonl")):
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
