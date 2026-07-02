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
