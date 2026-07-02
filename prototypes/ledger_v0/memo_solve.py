# This file is part of Bertini 2 (prototypes/ledger_v0 -- experimental, unshipped).
# GPL v3+; see the repository's licenses/ directory.

"""ensure_solved: solve() as "ensure this ask is answered" -- memoized against a Ledger.

- complete run recorded -> return it (no tracking);
- partial run -> adopt its recorded homotopy instance and finish the missing paths;
- nothing -> create the run (fresh randomness), record everything as it completes.

v0 stopgap: the homotopy instance (gamma, start coefficients) is persisted as a pickled
blob and re-adopted on resume -- the manifest mechanism, standing in until seed-rooted
randomness derivation exists (arc rung 2).
"""

import hashlib
import json
import pickle

import numpy as np

import bertini as pb
from bertini import multiprec as mp
from bertini.system import start_system as ss
from bertini import nag_algorithm as nag

from ledger import Ledger, SCHEMA


class SimulatedCrash(RuntimeError):
    """Raised by ensure_solved(crash_after=N) after N paths, standing in for a walltime kill."""


def _ask_digest(target, config: dict) -> dict:
    """The ask: WHAT was requested -- target identity + config.  (Seed joins in rung 2.)"""
    return {
        "target": target.content_digest(),
        "config": dict(sorted(config.items())),
    }


def _run_id(ask: dict) -> str:
    return hashlib.sha256(json.dumps(ask, sort_keys=True).encode()).hexdigest()[:16]


def _encode_point(vec) -> list:
    """Exact-ish text encoding of an endpoint: [re_str, im_str] per coordinate, at the
    solution's own precision (str of an mpfr carries its digits)."""
    out = []
    for z in vec:
        out.append([str(z.real), str(z.imag)])
    return out


def _decode_point(coords: list) -> np.ndarray:
    return np.array([mp.Complex(re, im) for re, im in coords])


class LedgerSolveResult:
    """What ensure_solved returns: solutions plus a small accounting of reuse."""

    def __init__(self, run_id, solutions, statuses, num_reused, num_computed):
        self.run_id = run_id
        self.solutions = solutions      # index -> np.ndarray (successful paths only)
        self.statuses = statuses        # index -> "success" | "failed"
        self.num_reused = num_reused
        self.num_computed = num_computed

    def all_solutions(self):
        return [self.solutions[i] for i in sorted(self.solutions)]


DEFAULT_CONFIG = {"precision": "adaptive", "endgame": "cauchy"}


def ensure_solved(target, ledger: Ledger, config=None, crash_after=None):
    """Ensure `target`'s total-degree zero-dim solve is answered in `ledger`."""
    config = dict(DEFAULT_CONFIG if config is None else config)
    ask = _ask_digest(target, config)
    run = ledger.find_run(ask)

    if run is None:
        run = _create_run(target, ledger, ask)

    run_id = run["run"]
    done = ledger.completed_paths(run_id)
    missing = [i for i in range(run["num_paths"]) if i not in done]

    num_computed = 0
    if missing:
        homotopy = pickle.loads(ledger.get_object(run["homotopy_blob"]))
        with ledger.open_journal(run_id) as journal:
            for i in missing:
                record = _track_one(target, homotopy, _decode_point(run["start_points"][i]), i, config)
                record["run"] = run_id
                record["start"] = {"kind": "start_label", "index": i}   # provenance bottoms out here
                journal.append(record)
                done[i] = record
                num_computed += 1
                if crash_after is not None and num_computed >= crash_after:
                    raise SimulatedCrash(
                        "simulated walltime kill after %d paths (run %s)" % (num_computed, run_id))

    solutions = {i: _decode_point(rec["endpoint"])
                 for i, rec in done.items() if rec["status"] == "success"}
    statuses = {i: rec["status"] for i, rec in done.items()}
    return LedgerSolveResult(run_id, solutions, statuses,
                             num_reused=len(done) - num_computed, num_computed=num_computed)


def _create_run(target, ledger: Ledger, ask: dict) -> dict:
    """Draw the instance (start system + homotopy, with their fresh randomness), persist
    it to the object store, and append the run header."""
    start = ss.TotalDegreeBinomial(target)
    homotopy = pb.system.make_homotopy(target, start, "t", None)  # None -> random gamma

    target_object = ledger.put_object(target.to_classic_input(), object_id=target.content_digest())
    homotopy_blob = ledger.put_object(pickle.dumps(homotopy))

    # the start points ARE provenance data (the canonical start labels' values), so they
    # live in the run header as text -- no need to persist the start-system object itself
    num_paths = start.num_start_points()
    start_points = [_encode_point(start.start_point_mp(i)) for i in range(num_paths)]

    run = {
        "kind": "run",
        "schema": SCHEMA,
        "run": _run_id(ask),
        "ask": ask,
        "target_object": target_object,
        "homotopy_blob": homotopy_blob,
        "num_paths": num_paths,
        "start_points": start_points,
    }
    with ledger.open_journal(run["run"]) as journal:
        journal.append(run)
    return run


def ensure_continued(target, generic, ledger: Ledger, config=None, crash_after=None):
    """The chain link: continue a previously-solved `generic` system's endpoints to
    `target` via the (deterministic, gamma=1) coefficient parameter homotopy.

    The generic's solve must already be on record (ensure_solved it first); its recorded
    endpoints -- read back from the records, not from memory -- become this run's start
    points, and each track record carries a `start` reference into the generic run:
    a real two-link provenance chain, walkable back to the total-degree start labels.

    No randomness enters here (gamma = 1), so the continuation homotopy is REBUILT from
    target+generic rather than persisted: chained runs need no blob at all.
    """
    config = dict(DEFAULT_CONFIG if config is None else config)
    ask = {
        "op": "continue",
        "target": target.content_digest(),
        "generic": generic.content_digest(),
        "config": dict(sorted(config.items())),
    }
    run = ledger.find_run(ask)

    if run is None:
        # any run that SOLVED the generic will do -- a base solve or itself a
        # continuation (chains nest: sample <- midpoint slice <- witness <- start)
        generic_run = _find_run_solving(ledger, generic.content_digest())
        if generic_run is None:
            raise ValueError("ensure_continued: the generic system has no recorded solve; "
                             "ensure_solved(generic, ledger) first")
        generic_done = ledger.completed_paths(generic_run["run"])
        start_indices = sorted(i for i, rec in generic_done.items() if rec["status"] == "success")

        target_object = ledger.put_object(target.to_classic_input(),
                                          object_id=target.content_digest())
        run = {
            "kind": "run",
            "schema": SCHEMA,
            "run": _run_id(ask),
            "ask": ask,
            "op": "continue",
            "target_object": target_object,
            "start_run": generic_run["run"],
            "start_indices": start_indices,
            "num_paths": len(start_indices),
        }
        with ledger.open_journal(run["run"]) as journal:
            journal.append(run)

    run_id = run["run"]
    done = ledger.completed_paths(run_id)
    missing = [i for i in range(run["num_paths"]) if i not in done]

    num_computed = 0
    if missing:
        homotopy = nag.coefficient_parameter_homotopy(target, generic)  # deterministic
        start_records = ledger.completed_paths(run["start_run"])
        with ledger.open_journal(run_id) as journal:
            for i in missing:
                generic_index = run["start_indices"][i]
                start_point = _decode_point(start_records[generic_index]["endpoint"])
                record = _track_one(target, homotopy, start_point, i, config)
                record["run"] = run_id
                record["start"] = {"kind": "point_ref",
                                   "run": run["start_run"], "index": generic_index}
                journal.append(record)
                done[i] = record
                num_computed += 1
                if crash_after is not None and num_computed >= crash_after:
                    raise SimulatedCrash(
                        "simulated walltime kill after %d paths (run %s)" % (num_computed, run_id))

    solutions = {i: _decode_point(rec["endpoint"])
                 for i, rec in done.items() if rec["status"] == "success"}
    statuses = {i: rec["status"] for i, rec in done.items()}
    return LedgerSolveResult(run_id, solutions, statuses,
                             num_reused=len(done) - num_computed, num_computed=num_computed)


def _find_run_solving(ledger: Ledger, target_digest: str):
    """The most recent run (of any op) whose ask.target is this system."""
    found = None
    for rec in ledger.scan():
        if rec.get("kind") == "run" and rec.get("ask", {}).get("target") == target_digest:
            found = rec
    return found


def annotate(ledger: Ledger, run_id: str, index: int, key: str, value):
    """Attach metadata (e.g. a projection value) to a recorded point: one more
    append-only line, keyed by point id -- queryable with jq/pandas, mergeable by
    file append."""
    with ledger.open_journal(run_id) as journal:
        journal.append({"kind": "annotation", "point": {"run": run_id, "index": index},
                        "key": key, "value": value})


def annotations_for(ledger: Ledger, run_id: str, index: int) -> dict:
    """All annotations recorded against one point, as {key: value}."""
    out = {}
    for rec in ledger.scan():
        if rec.get("kind") == "annotation" and rec.get("point") == {"run": run_id, "index": index}:
            out[rec["key"]] = rec["value"]
    return out


def provenance_chain(ledger: Ledger, run_id: str, index: int) -> list:
    """Walk one endpoint's ancestry back to the beginning: the arc's 'all the way to the
    start' promise, executed against nothing but the records."""
    chain = []
    while True:
        rec = ledger.completed_paths(run_id).get(index)
        if rec is None:
            raise KeyError("no track record for (%s, %d)" % (run_id, index))
        chain.append({"run": run_id, "index": index, "status": rec["status"]})
        start = rec.get("start", {})
        if start.get("kind") == "point_ref":
            run_id, index = start["run"], start["index"]   # ascend one link
        else:
            chain.append({"run": run_id, "start_label": start.get("index", index)})
            return chain


def _track_one(target, homotopy, point, index: int, config: dict) -> dict:
    """Track a single start point; one path per solver call keeps index -> endpoint exact."""
    solver = nag.HomotopySolver(homotopy, [point], target,
                                precision=config["precision"], endgame=config["endgame"])
    solver.solve()
    endpoints = solver.all_solutions()
    if len(endpoints) == 1:
        return {"kind": "track", "index": index, "status": "success",
                "endpoint": _encode_point(endpoints[0])}
    return {"kind": "track", "index": index, "status": "failed", "endpoint": None}
