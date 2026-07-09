# This file is part of Bertini 2.
#
# python/bertini/records.py is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Bertini 2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
# See the GNU General Public License for more details.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.

"""The casual face of the structured output directory: ``solve``, ``save``, ``load``.

Every solve writes durable, plain-text records of what was computed (the directory
explains itself: see its README.txt), and consults them before computing -- so rerunning
a script is always safe: instant if done, resuming if crashed, fresh if new.  The
records directory is ambient (``BERTINI_RECORDS_DIR``, else ``./bertini_output``);
nobody is required to name it.

    import bertini as pb
    sols = pb.solve(my_system, seed=42)   # records + resumes automatically
    pb.save(sols)                         # or pb.save("my favorites", sols)
    pb.load("my favorites")               # back, in any later session

Reading the records never requires bertini: JSON lines throughout -- ``history/``
(what was asked, when), ``results/`` (one file per run: the computed paths and their
metadata), ``definitions/`` (the exact inputs, content-addressed).  Power users keep
the full solver-object API; these three verbs are sugar over it.
"""

import json as _json
import os as _os
import time as _time
from pathlib import Path as _Path

import numpy as _np

__all__ = ['solve', 'save', 'load', 'records_dir', 'Solution', 'SolveResult']


# --- the ambient directory -----------------------------------------------------------

_ambient = None


_recording_enabled = True


def recording(on=None):
    """The one-line on/off switch for the records: ``bertini.recording(False)``.

    Recording is ON BY DEFAULT for every solver in the process -- ``solve``,
    ``ZeroDimSolver``, ``HomotopySolver`` -- into the ambient directory (see
    ``records_dir``).  With recording off, solves run bare: no directory is written
    or consulted (so no resuming either), and the ambient ``BERTINI_RECORDS_DIR``
    attach is suppressed for this process.  ``recording(True)`` turns it back on.
    Call with no argument to ask the current state.

    For the command line, the same switch is the environment:
    ``BERTINI_RECORDS_DIR=none`` runs ``bertini2`` without records (the empty string
    also works on POSIX; Windows deletes empty-valued variables, so ``none`` is the
    portable spelling).
    """
    global _recording_enabled
    if on is not None:
        _recording_enabled = bool(on)
        # the C++ side reads the environment for its ambient attach; keep it in step
        # ('none' rather than '' -- Windows deletes empty-valued variables, and with
        # records on by default a deleted variable would mean ON)
        if not _recording_enabled:
            _os.environ['BERTINI_RECORDS_DIR'] = 'none'
        elif _ambient is not None:
            _os.environ['BERTINI_RECORDS_DIR'] = _ambient
        elif _os.environ.get('BERTINI_RECORDS_DIR') in ('', 'none'):
            del _os.environ['BERTINI_RECORDS_DIR']
    return _recording_enabled


def records_dir(path=None):
    """Get (or set, by passing a path) the ambient records directory for this process.

    Resolution when unset: the ``BERTINI_RECORDS_DIR`` environment variable, else
    ``./bertini_output``.  Created on demand.  Returns the resolved path as a string.

    Ambient recording is ON BY DEFAULT for every solver in the process --
    ``ZeroDimSolver`` and ``HomotopySolver`` record here just like ``solve``, with no
    code at all.  Setting a path chooses WHERE: it is exported to
    ``BERTINI_RECORDS_DIR``, which the solver classes read for their ambient attach.
    ``recording(False)`` is the off switch: while recording is off, nothing is
    exported and the off sentinel survives.
    """
    global _ambient
    if path is not None:
        _ambient = str(path)
        # the solver classes' ambient attach reads the environment; keep it in step
        # (unless recording is switched off, whose off sentinel must survive)
        if _recording_enabled:
            _os.environ['BERTINI_RECORDS_DIR'] = _ambient
    if _ambient is None:
        env = _os.environ.get('BERTINI_RECORDS_DIR')
        # '' / 'none' are the recording-off sentinels, not directory names
        _ambient = env if env and env != 'none' else 'bertini_output'
    return _ambient


def _directory(path=None):
    """The bound OutputDirectory writer (the single C++ implementation) at `path`."""
    from bertini._pybertini.records import OutputDirectory
    return OutputDirectory(str(path if path is not None else records_dir()))


# --- solutions are points that remember ----------------------------------------------

class Solution(_np.ndarray):
    """A solution point: coordinates that remember where they came from.

    Behaves exactly like the numpy array you expect (index it, print it, feed it to a
    solver), while carrying ``.provenance`` (``{'run': ..., 'index': ...}`` -- the
    recorded path that produced it) and ``.annotations`` invisibly.  Arithmetic
    produces plain derived points: a computed combination is a new thing, and its
    provenance is honestly absent.
    """

    def __new__(cls, coordinates, provenance=None, annotations=None):
        obj = _np.asarray(coordinates).view(cls)
        obj.provenance = provenance
        obj.annotations = dict(annotations or {})
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        # views keep memory; arithmetic results get fresh (empty) provenance
        self.provenance = getattr(obj, 'provenance', None)
        self.annotations = dict(getattr(obj, 'annotations', {}) or {})

    def __array_wrap__(self, out_arr, context=None, return_scalar=False):
        # a derived point is a NEW point: provenance does not survive arithmetic
        result = super().__array_wrap__(out_arr, context, return_scalar)
        if isinstance(result, Solution) and context is not None:
            result.provenance = None
            result.annotations = {}
        return result

    # numpy's ndarray .real/.imag are hardwired to its built-in complex types: on a
    # multiprecision-complex array the base attributes return silently WRONG values
    # (.real gives the complex values, .imag gives zeros), with no hook for user
    # dtypes.  A python subclass property CAN shadow the C-level attribute, so
    # solution points -- the mp-complex arrays users actually hold -- are correct.
    # Plain ndarrays are covered by the guarded np.real/np.imag/np.angle (see
    # bertini._numpy_guard) and by bertini.real/imag.

    @property
    def real(self):
        """The real parts -- correct also for the multiprecision complex dtype."""
        from bertini.multiprec import complex_mp
        if self.dtype == _np.dtype(complex_mp):
            from bertini import _numpy_helpers as _nh
            return _nh.real(_np.asarray(self))
        return _np.ndarray.real.__get__(self)

    @property
    def imag(self):
        """The imaginary parts -- correct also for the multiprecision complex dtype."""
        from bertini.multiprec import complex_mp
        if self.dtype == _np.dtype(complex_mp):
            from bertini import _numpy_helpers as _nh
            return _nh.imag(_np.asarray(self))
        return _np.ndarray.imag.__get__(self)


class SolveResult:
    """What ``solve`` returns: the solutions plus a claim ticket on the recorded run.

    Forgetting to capture it loses nothing -- the records hold the truth; another
    ``solve`` of the same ask re-mints an equivalent result (recalled, not recomputed).
    """

    def __init__(self, solutions, run_id, directory, num_recalled, solver):
        self.solutions = solutions          #: list[Solution]: the finite solutions, user coordinates
        self.run_id = run_id                #: str: the recorded run's id ({run, index} is a point reference)
        self.directory = directory          #: str: the records directory this run lives in
        self.num_recalled = num_recalled    #: int: paths taken from the records instead of computed
        self._solver = solver               # kept alive: the power-user escape hatch

    def __len__(self):
        return len(self.solutions)

    def __iter__(self):
        return iter(self.solutions)

    def __getitem__(self, k):
        return self.solutions[k]

    def __repr__(self):
        return ("SolveResult(%d solutions, run %s, %d recalled, records at %s)"
                % (len(self.solutions), self.run_id, self.num_recalled, self.directory))

    @property
    def solver(self):
        """The underlying solver object (all_solutions, solution_metadata, ...)."""
        return self._solver


# --- chained solves -------------------------------------------------------------------

def _point_reference(point, position, pending_givens):
    """The provenance reference for one start point: a point_ref when the point carries
    ``.provenance`` {run, index} (a Solution from a prior solve -- a CHAIN), else a
    slot in the run's given (external data -- provenance bottoms out honestly)."""
    prov = getattr(point, 'provenance', None)
    if isinstance(prov, dict) and 'run' in prov and 'index' in prov:
        return {'kind': 'point_ref', 'run': str(prov['run']), 'index': int(prov['index'])}
    pending_givens.append((position, point))
    return None    # patched once the given is archived


def _archive_given_points(pending, directory_writer):
    """Archive externally supplied start points as ONE given definition (coordinates as
    full-precision strings, readable without bertini) plus its given record; returns the
    definition id."""
    rows = []
    for _, point in pending:
        coords = []
        for v in _np.atleast_1d(_np.asarray(point, dtype=object)).ravel():
            c = complex(v)
            coords.append([repr(c.real), repr(c.imag)])
        rows.append(coords)
    content = _json.dumps({'kind': 'start_points', 'points': rows}, indent=1)
    given_id = directory_writer.put_definition(content, 'givens', label='start_points')
    directory_writer.append(_json.dumps({
        'kind': 'given', 'source': given_id, 'role': 'start_points',
        'when': _time.strftime('%Y-%m-%d %H:%M')}))
    return given_id


def _coerce_start_point(point, precision):
    """Coerce one start point to the coordinate type the solver's converter expects:
    multiprecision complex for 'adaptive'/'multiple' (a raw float/complex array would
    fail eigenpy conversion), plain complex for 'double'.  Points already holding mp
    coordinates (a prior solve's solutions) pass through untouched."""
    arr = _np.atleast_1d(_np.asarray(point))
    if precision == 'double':
        return arr.astype(complex)
    from bertini.multiprec import complex_mp
    mp_dtype = _np.dtype(complex_mp)   # the eigenpy-registered numpy dtype
    if arr.dtype == mp_dtype:
        return arr
    return _np.array([v if isinstance(v, complex_mp)
                      else complex_mp(repr(complex(v).real), repr(complex(v).imag))
                      for v in arr.ravel()], dtype=mp_dtype)


def _chained_solver(system, homotopy, start, where, *, precision, endgame):
    """Build the HomotopySolver for a chained solve, plus the per-path provenance refs
    and the start-data identity that joins the ask."""
    import hashlib
    from bertini import nag_algorithm as _nag

    points = list(getattr(start, 'solutions', start))
    if not points:
        raise ValueError('solve: start= supplied no points')

    pending_givens = []
    refs = [_point_reference(pt, k, pending_givens) for k, pt in enumerate(points)]
    if pending_givens:
        given_id = _archive_given_points(pending_givens, _directory(where))
        for slot, (position, _) in enumerate(pending_givens):
            refs[position] = {'kind': 'given_ref', 'given': given_id, 'index': slot}

    # the identity of the start data joins the ask: the same homotopy from different
    # start points is a different computation
    identity = hashlib.sha256(
        _json.dumps(refs, sort_keys=True).encode()).hexdigest()

    solver = _nag.HomotopySolver(homotopy,
                                 [_coerce_start_point(pt, precision) for pt in points],
                                 system, precision=precision, endgame=endgame)
    return solver, refs, identity


# --- the three verbs ------------------------------------------------------------------

def solve(system, seed=None, directory=None, precision='adaptive', endgame='cauchy',
          homotopy=None, start=None):
    """Solve a polynomial system, recording and resuming automatically.

    Ensure-answered semantics: the solve consults the ambient records directory first;
    paths already recorded for this exact ask (system + settings + seed) are taken from
    the records, and only the rest are computed.  Rerun a crashed script and it
    finishes; rerun a finished one and it is instant.

    Parameters
    ----------
    system : System
        The target system (no path variable).
    seed : int, optional
        The reproducibility seed.  ``solve(sys, seed=42)`` means the same homotopy --
        the same gamma, start points, and patch -- forever, on every machine.  Omitted:
        an EFFECTIVE seed is derived for this solve and the solve runs under it, so
        the seed recorded in the run's ask reproduces that run standalone -- a
        mid-session solve does not silently depend on the session's earlier draw
        history.  (Consecutive seedless solves get distinct seeds, chained
        deterministically from ``set_random_seed``'s master when one was set.)
    directory : str, optional
        Records directory override; default is ambient (see ``records_dir``).
    homotopy : System, optional
        A homotopy you built (e.g. :func:`bertini.nag_algorithm.blend_homotopy`), for a
        CHAINED solve: its paths run from your ``start`` points at t=1 to ``system``'s
        solutions at t=0.  Requires ``start``.
    start : SolveResult or iterable of points, optional
        Where the paths start.  A prior :class:`SolveResult` (or its solutions) chains
        with full provenance -- the records link every new endpoint back through the
        prior run, all the way to the beginning.  Raw points (arrays) are archived as a
        *given*: provenance bottoms out honestly at data you supplied.
    precision, endgame : str
        Passed through to :func:`bertini.nag_algorithm.ZeroDimSolver` (``mptype`` /
        ``endgame``).

    Returns
    -------
    SolveResult
        The finite solutions (as :class:`Solution` points that remember their run) plus
        the run id -- a claim ticket, safe to drop.
    """
    from bertini import nag_algorithm as _nag
    from bertini.random import set_random_seed as _set_seed, derive_solve_seed as _derive_seed

    if (homotopy is None) != (start is None):
        raise ValueError("solve: homotopy= and start= go together (a chained solve "
                         "needs both the homotopy and where its paths start)")
    if seed is None:
        # derive this solve's own effective seed: the run's recorded seed must
        # reproduce the run STANDALONE, never depend on the session's earlier draw
        # history.  (A chained solve's homotopy was built before this line -- its
        # exact coefficients are archived in definitions/ by the run header, so
        # nothing is lost there either.)
        seed = _derive_seed()
    _set_seed(seed)

    where = str(directory if directory is not None else records_dir())
    if homotopy is not None:
        zd, refs, identity = _chained_solver(system, homotopy, start, where,
                                             precision=precision, endgame=endgame)
        if _recording_enabled:
            zd.record_to(where)
            zd.set_recorded_start_provenance(_json.dumps(refs), identity)
    else:
        zd = _nag.ZeroDimSolver(system, mptype=precision, endgame=endgame)
        if _recording_enabled:
            zd.record_to(where)
    zd.solve()

    run_id = zd.records_run_id()
    # finite solutions, carrying their TRUE path indices as provenance ({run, index}
    # is exactly how the records reference points); with recording off there is no
    # run to reference, so provenance is honestly absent
    all_sols = zd.all_solutions()
    solutions = [Solution(all_sols[int(m.path_index)],
                          provenance=({'run': run_id, 'index': int(m.path_index)}
                                      if run_id else None))
                 for m in zd.solution_metadata() if m.is_finite]

    result = SolveResult(solutions, run_id, where if run_id else None,
                         int(zd.num_paths_recalled()), zd)
    if run_id:
        # top-level solves auto-declare their deliverable: "what were my solutions?"
        save("solutions [run %s]" % run_id, result,
             description="finite solutions, auto-declared by bertini.solve",
             directory=where)
    return result


def save(*args, description='', directory=None):
    """Save almost anything under a name: ``save(thing)`` or ``save(name, thing)``.

    A :class:`SolveResult` (or anything with ``.run_id`` and ``.solutions``) is declared
    as results with full provenance; any JSON-able value (dict, list, number, string) is
    recorded inline.  The declaration lands in ``history/`` (the points themselves
    live in ``results/``, referred to by {run, index}).  A nameless ``save(thing)``
    is auto-named by timestamp; re-saving a name replaces it (newest wins).
    """
    if len(args) == 1:
        name, thing = 'saved %s' % _time.strftime('%Y-%m-%d %H:%M:%S'), args[0]
    elif len(args) == 2:
        name, thing = args
    else:
        raise TypeError('save(thing) or save(name, thing)')

    record = {'kind': 'result', 'name': str(name), 'description': description,
              'when': _time.strftime('%Y-%m-%d %H:%M')}
    if hasattr(thing, 'run_id') and hasattr(thing, 'solutions'):
        record['points'] = [{'run': thing.run_id, 'index': (s.provenance or {}).get('index', i)}
                            for i, s in enumerate(thing.solutions)]
    else:
        record['points'] = []
        record['value'] = thing   # any JSON-able thing

    out = _directory(directory)
    out.append(_json.dumps(record))
    return name


def _scan_history(directory=None):
    """Every record in the directory's history, in journal order (plain-json read)."""
    root = _Path(directory if directory is not None else records_dir()) / 'history'
    records = []
    if root.is_dir():
        for journal in sorted(root.glob('*.jsonl')):
            for line in journal.read_text().splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    records.append(_json.loads(line))
                except ValueError:
                    continue    # a torn final line from a killed writer
    return records


def _scan_results(directory=None, run=None):
    """Every record in the directory's results files (plain-json read).

    ``run=`` reads just that run's file (``results/<2 hex>/<run>.jsonl``).  Header
    lines (kind ``results_header``) are included; callers filter by kind.
    """
    root = _Path(directory if directory is not None else records_dir()) / 'results'
    if run is not None:
        run = str(run)
        path = root / run[:2] / (run + '.jsonl')
        files = [path] if path.exists() else []
    else:
        files = sorted(root.glob('*/*.jsonl')) if root.is_dir() else []
    records = []
    for path in files:
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                records.append(_json.loads(line))
            except ValueError:
                continue    # a torn final line from a killed writer
    return records


def solutions_of(run, directory=None, status='success'):
    """The recorded endpoints of a run, as :class:`Solution` points with provenance.

    Works on ANY run in the directory -- including runs written by the command-line
    ``bertini2`` or on another machine: this reads the plain records, no solver object
    needed.  The points are in USER coordinates, exact to the recorded precision, and
    carry ``{'run', 'index'}`` provenance, so they chain directly:
    ``solve(B, homotopy=H, start=solutions_of(run_id))``.

    ``status`` filters by recorded verdict (``'success'`` default; pass ``None`` for
    every tracked path -- audits want all three populations).
    """
    from bertini.multiprec import complex_mp
    tracks = {}
    for rec in _scan_results(directory, run=run):
        if rec.get('kind') == 'path' and 'index' in rec:
            tracks[int(rec['index'])] = rec      # newest wins per index
    points = []
    for index in sorted(tracks):
        rec = tracks[index]
        if status is not None and rec.get('status') != status:
            continue
        endpoint = rec.get('endpoint_user') or rec.get('endpoint')
        if not endpoint:
            continue
        coords = _np.array([complex_mp(str(c[0]), str(c[1])) for c in endpoint],
                           dtype=_np.dtype(complex_mp))
        points.append(Solution(coords, provenance={'run': str(run), 'index': index}))
    return points


def provenance(point, directory=None):
    """Walk a point's provenance back to the beginning: the chain of ``{'run','index'}``
    hops, ending at a canonical start label or a given (user-supplied data).

    ``point`` is a :class:`Solution` (or anything with ``.provenance``), or an explicit
    ``{'run': ..., 'index': ...}`` dict.  Reads the plain records -- runs written by
    the CLI and by Python walk the same.
    """
    ref = getattr(point, 'provenance', None) or point
    if not isinstance(ref, dict) or 'run' not in ref or 'index' not in ref:
        raise ValueError("provenance needs a point with {'run', 'index'}")
    tracks = {}
    for rec in _scan_results(directory):
        if rec.get('kind') == 'path' and 'run' in rec and 'index' in rec:
            tracks[(str(rec['run']), int(rec['index']))] = rec
    hops = []
    current = {'run': str(ref['run']), 'index': int(ref['index'])}
    while True:
        hops.append(current)
        rec = tracks.get((current['run'], current['index']))
        if rec is None:
            hops.append({'kind': 'unrecorded'})
            return hops
        start = rec.get('start', {})
        if start.get('kind') == 'point_ref':
            current = {'run': str(start['run']), 'index': int(start['index'])}
            continue
        hops.append(start)          # start_label or given_ref: the beginning
        return hops


# --- navigating the records: tables and the provenance graph ---------------------------
#
# Scale is a design constraint here: a directory may hold MILLIONS of tracked paths.
# The tables stream plain records into pandas (fine at that scale); coordinates stay
# OUT of the tables unless asked (they are the bulk of the bytes); and the chain plot
# aggregates to run level rather than drawing a node per path once a run is large.

def runs(directory=None):
    """One row per recorded run, as a :class:`pandas.DataFrame`.

    Columns: ``run``, ``when``, ``op``, ``num_paths``, ``seed``, ``tracker``,
    ``endgame``, ``target_digest``, ``recalls`` (how many later sessions answered
    this ask from the store -- 0 means computed once, never re-asked),
    ``producer_version``, ``producer_commit``.
    Reads the plain records -- any directory, any producer (CLI or Python).
    """
    import pandas as pd
    history = _scan_history(directory)
    recall_counts = {}
    for rec in history:
        if rec.get('kind') == 'recall':
            run_id = rec.get('run')
            recall_counts[run_id] = recall_counts.get(run_id, 0) + 1
    rows = []
    for rec in history:
        if rec.get('kind') != 'run':
            continue
        ask = rec.get('ask', {})
        producer = rec.get('producer', {})
        rows.append({
            'run': rec.get('run'),
            'when': rec.get('when'),
            'op': rec.get('op'),
            'num_paths': rec.get('num_paths'),
            'seed': ask.get('seed'),
            'tracker': ask.get('tracker'),
            'endgame': ask.get('endgame'),
            'target_digest': rec.get('target_digest'),
            'recalls': recall_counts.get(rec.get('run'), 0),
            'producer_version': producer.get('version'),
            'producer_commit': producer.get('commit'),
        })
    return pd.DataFrame(rows)


def tracks(run=None, directory=None, coordinates=False):
    """One row per tracked path (newest record per (run, index)), as a DataFrame.

    Columns: ``run``, ``index``, ``status``, ``outcome`` (the endgame verdict's name),
    ``start_kind``, ``start_run``, ``start_index``, ``cycle_num``,
    ``path_time_seconds``.  Pass ``run=`` to restrict to one run.

    ``coordinates=False`` (the default) keeps endpoints OUT of the table -- they are
    most of the bytes, and a million-path audit usually wants the statuses, not the
    numbers.  ``coordinates=True`` adds an ``endpoint_user`` column of complex tuples.
    """
    import pandas as pd
    newest = {}
    for rec in _scan_results(directory, run=run):
        if rec.get('kind') != 'path' or 'run' not in rec or 'index' not in rec:
            continue
        newest[(rec['run'], int(rec['index']))] = rec
    rows = []
    for (run_id, index), rec in sorted(newest.items()):
        start = rec.get('start', {})
        row = {
            'run': run_id,
            'index': index,
            'status': rec.get('status'),
            'outcome': rec.get('endgame_success_code_name'),
            'start_kind': start.get('kind'),
            'start_run': start.get('run'),
            'start_index': start.get('index'),
            'cycle_num': rec.get('cycle_num'),
            'path_time_seconds': (float(rec['path_time_seconds'])
                                  if 'path_time_seconds' in rec else None),
        }
        if coordinates:
            endpoint = rec.get('endpoint_user') or rec.get('endpoint') or []
            row['endpoint_user'] = tuple(complex(float(c[0]), float(c[1]))
                                         for c in endpoint)
        rows.append(row)
    return pd.DataFrame(rows)


def provenance_graph(directory=None, runs=None):
    """The provenance of every recorded point, as a :class:`networkx.DiGraph`.

    Nodes are points ``('run_id', index)`` plus origins ``('start_label', i)`` /
    ``('given', definition_id, i)``; each edge points FROM where a path started TO its
    endpoint (time flows along edges).  Node attributes: ``run``, ``index``,
    ``status``; edge attribute ``run`` (the run that tracked it).

    Pass ``runs=`` (an iterable of run ids) to restrict a large directory to the
    chains you care about -- a million-path directory makes a million-node graph,
    which networkx holds but no drawing survives (see :func:`plot_chain`, which
    aggregates instead).
    """
    import networkx as nx
    wanted = set(runs) if runs is not None else None
    graph = nx.DiGraph()
    for rec in _scan_results(directory):
        if rec.get('kind') != 'path' or 'run' not in rec or 'index' not in rec:
            continue
        run_id = str(rec['run'])
        if wanted is not None and run_id not in wanted:
            continue
        index = int(rec['index'])
        endpoint = (run_id, index)
        graph.add_node(endpoint, run=run_id, index=index, status=rec.get('status'))
        start = rec.get('start', {})
        kind = start.get('kind')
        if kind == 'point_ref':
            origin = (str(start['run']), int(start['index']))
        elif kind == 'given_ref':
            origin = ('given', str(start.get('given', '')), int(start.get('index', -1)))
        else:
            origin = ('start_label', run_id, int(start.get('index', index)))
        graph.add_edge(origin, endpoint, run=run_id)
    return graph


def plot_chain(directory=None, runs=None, ax=None, max_paths_drawn=200):
    """Draw the chain left to right: each run is a column, paths flow rightward from
    their starts to their endpoints.

    Small chains draw every path (green = success, orange = diverged, red = failed;
    origins are squares).  A run with more than ``max_paths_drawn`` paths is NOT drawn
    path-by-path -- the whole figure falls back to one node per run with edge widths
    showing how many paths flow between runs, so a million-path chain renders in
    milliseconds instead of crashing your session.

    Returns the matplotlib Axes.
    """
    import networkx as nx
    import matplotlib.pyplot as plt

    graph = provenance_graph(directory, runs=runs)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 5))

    # column per run, ordered by chain depth: origins at 0, then each run one right of
    # the deepest run it draws starts from
    run_edges = {}
    run_sizes = {}
    for origin, endpoint, data in graph.edges(data=True):
        run_id = data['run']
        run_sizes[run_id] = run_sizes.get(run_id, 0) + 1
        if origin in graph.nodes and len(origin) == 2:
            source = graph.nodes[origin].get('run')
            if source is not None and source != run_id:
                run_edges.setdefault(run_id, set()).add(source)
    depth = {}
    def run_depth(run_id, _seen=()):
        if run_id in depth:
            return depth[run_id]
        parents = run_edges.get(run_id, ())
        d = 1 + max((run_depth(p) for p in parents if p not in _seen), default=0)
        depth[run_id] = d
        return d
    for run_id in run_sizes:
        run_depth(run_id)

    if run_sizes and max(run_sizes.values()) > max_paths_drawn:
        # AGGREGATE: one node per run; edges weighted by path counts between runs
        agg = nx.DiGraph()
        for run_id, size in run_sizes.items():
            agg.add_node(run_id, subset=depth[run_id], size=size)
        counts = {}
        for origin, endpoint, data in graph.edges(data=True):
            source = graph.nodes.get(origin, {}).get('run')
            if source is not None and source != data['run']:
                counts[(source, data['run'])] = counts.get((source, data['run']), 0) + 1
        for (a, b), n in counts.items():
            agg.add_edge(a, b, weight=n)
        pos = nx.multipartite_layout(agg, subset_key='subset')
        widths = [1 + 4 * agg[a][b]['weight'] / max(counts.values()) for a, b in agg.edges]
        nx.draw_networkx(agg, pos, ax=ax, node_color='#88aadd', node_size=900,
                         width=widths, font_size=7,
                         labels={r: '%s\n(%d paths)' % (r[:8], run_sizes[r]) for r in agg})
        ax.set_title('chain of runs (aggregated: runs exceed max_paths_drawn)')
        import matplotlib.lines as mlines
        ax.legend(handles=[mlines.Line2D([], [], color='#777777', linewidth=3,
                                         label='paths flowing between runs (width = count)')],
                  loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    else:
        for node in graph.nodes:
            data = graph.nodes[node]
            graph.nodes[node]['subset'] = (depth[data['run']] if data.get('run') in depth
                                           else 0)
        pos = nx.multipartite_layout(graph, subset_key='subset')
        status_color = {'success': '#2a9d3a', 'diverged': '#e8930c', 'failed': '#d43a2f'}
        colors = [status_color.get(graph.nodes[n].get('status'), '#888888')
                  for n in graph.nodes]
        shapes = ['s' if graph.nodes[n].get('run') not in depth else 'o'
                  for n in graph.nodes]
        nx.draw_networkx_edges(graph, pos, ax=ax, arrows=True, edge_color='#777777')
        for shape in set(shapes):
            nodes = [n for n, sh in zip(graph.nodes, shapes) if sh == shape]
            nx.draw_networkx_nodes(graph, pos, nodelist=nodes, ax=ax, node_shape=shape,
                                   node_color=[colors[i] for i, n in enumerate(graph.nodes)
                                               if shapes[i] == shape], node_size=160)
        ax.set_title('paths through the chain')
        import matplotlib.lines as mlines
        legend_entries = [
            mlines.Line2D([], [], color='#aaaaaa', marker='s', linestyle='none',
                          markersize=8, label='start point (the beginning)'),
            mlines.Line2D([], [], color='#2a9d3a', marker='o', linestyle='none',
                          markersize=8, label='endpoint: success'),
            mlines.Line2D([], [], color='#e8930c', marker='o', linestyle='none',
                          markersize=8, label='endpoint: diverged (truncated)'),
            mlines.Line2D([], [], color='#d43a2f', marker='o', linestyle='none',
                          markersize=8, label='endpoint: failed'),
        ]
        ax.legend(handles=legend_entries, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                  fontsize=8, frameon=False)
    # the horizontal axis is the chain itself: solves left to right
    ax.set_axis_off()
    ax.annotate('chain depth: each column is one run, tracked left to right \u2192',
                xy=(0.5, 0.02), xycoords='axes fraction', ha='center', fontsize=9,
                color='#444444')
    return ax


def annotate(point, key, value, directory=None):
    """Attach metadata to a solution: ``annotate(sol, 'projection', 1.5)``.

    ``point`` is a :class:`Solution` (or anything with ``.provenance`` holding
    ``{'run', 'index'}``), or an explicit ``{'run': ..., 'index': ...}`` dict.  The
    annotation lands in the records beside the point it describes (readers such as
    :func:`load` merge it in); re-annotating the same key replaces it (newest wins).
    ``value`` is any JSON-able thing.
    """
    provenance = getattr(point, 'provenance', None) or point
    if not isinstance(provenance, dict) or 'run' not in provenance or 'index' not in provenance:
        raise ValueError("annotate needs a point with provenance {'run', 'index'} "
                         "(a Solution from solve(), or an explicit dict)")
    out = _directory(directory)
    out.annotate(str(provenance['run']), int(provenance['index']), str(key),
                 _json.dumps(value))
    if hasattr(point, 'annotations'):
        point.annotations[str(key)] = value


def load(name=None, directory=None):
    """Load saved results by name -- the other half of :func:`save`.

    ``load("my favorites")`` returns that result (its points, annotations, provenance,
    or its inline value); ``load()`` returns the whole dict of everything saved, by
    name.  Reads the plain records -- declarations and annotations from ``history/``,
    the referenced points from ``results/`` -- so this works in any later session, on
    any producer's directory, and the same files are readable without bertini at all.
    """
    runs_by_id = {}
    declared = {}
    annotations = {}
    for rec in _scan_history(directory):
        kind = rec.get('kind')
        if kind == 'run':
            runs_by_id[str(rec.get('run'))] = rec
        elif kind == 'result':
            declared[rec.get('name')] = rec          # newest declaration wins
        elif kind == 'annotation':
            pt = rec.get('point', {})
            where = (str(pt.get('run')), int(pt.get('index', -1)))
            annotations.setdefault(where, {})[rec.get('key')] = rec.get('value')

    paths = {}          # (run, index) -> newest payload record, loaded per run on demand
    loaded_runs = set()

    def _payload(run_id, index):
        if run_id not in loaded_runs:
            loaded_runs.add(run_id)
            for rec in _scan_results(directory, run=run_id):
                if rec.get('kind') == 'path' and 'index' in rec:
                    paths[(run_id, int(rec['index']))] = rec
        return paths.get((run_id, index))

    def _entry(rec):
        entry = {'declared': rec.get('when'), 'description': rec.get('description', '')}
        if 'value' in rec:
            entry['value'] = rec['value']
        points = []
        for ref in rec.get('points', []):
            run_id, index = str(ref.get('run')), int(ref.get('index', -1))
            point = {'provenance': ref}
            payload = _payload(run_id, index)
            if payload is None:
                point['status'] = 'missing'
                points.append(point)
                continue
            # the recorded verdict, verbatim: success / diverged / failed -- a
            # diverged (truncated-near-infinity) path is an answer, not a failure
            point['status'] = payload.get('status', 'success')
            if 'endgame_success_code_name' in payload:
                point['outcome'] = payload['endgame_success_code_name']
            # prefer the endpoint in USER coordinates (with the user's variable
            # names, from the run header); fall back to the internal endpoint
            endpoint = payload.get('endpoint_user') or payload.get('endpoint') or []
            names_key = 'variables_user' if 'endpoint_user' in payload else 'variables'
            var_names = (runs_by_id.get(run_id) or {}).get(names_key, [])
            if len(var_names) != len(endpoint):
                var_names = []
            point['coordinates'] = {
                (var_names[k] if k < len(var_names) else 'coordinate_%d' % k): c
                for k, c in enumerate(endpoint)}
            point['annotations'] = annotations.get((run_id, index), {})
            points.append(point)
        if points:
            entry['points'] = points
        return entry

    if name is not None:
        rec = declared.get(name)
        return None if rec is None else _entry(rec)
    return {n: _entry(r) for n, r in declared.items()}
