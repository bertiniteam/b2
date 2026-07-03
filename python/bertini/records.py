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

Reading the records never requires bertini: they are JSON lines (``history/``) and a
``results.json`` one ``json.load`` away.  Power users keep the full solver-object API;
these three verbs are sugar over it.
"""

import json as _json
import os as _os
import time as _time
from pathlib import Path as _Path

import numpy as _np

__all__ = ['solve', 'save', 'load', 'records_dir', 'Solution', 'SolveResult']


# --- the ambient directory -----------------------------------------------------------

_ambient = None


def records_dir(path=None):
    """Get (or set, by passing a path) the ambient records directory for this process.

    Resolution when unset: the ``BERTINI_RECORDS_DIR`` environment variable, else
    ``./bertini_output``.  Created on demand.  Returns the resolved path as a string.
    """
    global _ambient
    if path is not None:
        _ambient = str(path)
    if _ambient is None:
        _ambient = _os.environ.get('BERTINI_RECORDS_DIR', 'bertini_output')
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


class SolveResult:
    """What ``solve`` returns: the solutions plus a claim ticket on the recorded run.

    Forgetting to capture it loses nothing -- the records hold the truth; another
    ``solve`` of the same ask re-mints an equivalent result (hydrated, not recomputed).
    """

    def __init__(self, solutions, run_id, directory, num_hydrated, solver):
        self.solutions = solutions          #: list[Solution]: the finite solutions, user coordinates
        self.run_id = run_id                #: str: the recorded run's id ({run, index} is a point reference)
        self.directory = directory          #: str: the records directory this run lives in
        self.num_hydrated = num_hydrated    #: int: paths taken from the records instead of computed
        self._solver = solver               # kept alive: the power-user escape hatch

    def __len__(self):
        return len(self.solutions)

    def __iter__(self):
        return iter(self.solutions)

    def __getitem__(self, k):
        return self.solutions[k]

    def __repr__(self):
        return ("SolveResult(%d solutions, run %s, %d hydrated, records at %s)"
                % (len(self.solutions), self.run_id, self.num_hydrated, self.directory))

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
    given_id = directory_writer.put_definition(content)
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
        a fresh seed is drawn (and recorded in the run's ask).
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
    from bertini.random import set_random_seed as _set_seed

    if (homotopy is None) != (start is None):
        raise ValueError("solve: homotopy= and start= go together (a chained solve "
                         "needs both the homotopy and where its paths start)")
    if seed is not None:
        _set_seed(seed)

    where = str(directory if directory is not None else records_dir())
    if homotopy is not None:
        zd, refs, identity = _chained_solver(system, homotopy, start, where,
                                             precision=precision, endgame=endgame)
        zd.record_to(where)
        zd.set_recorded_start_provenance(_json.dumps(refs), identity)
    else:
        zd = _nag.ZeroDimSolver(system, mptype=precision, endgame=endgame)
        zd.record_to(where)
    zd.solve()
    zd.refresh_results()

    run_id = zd.records_run_id()
    # finite solutions, carrying their TRUE path indices as provenance ({run, index}
    # is exactly how the records reference points)
    all_sols = zd.all_solutions()
    solutions = [Solution(all_sols[int(m.path_index)],
                          provenance={'run': run_id, 'index': int(m.path_index)})
                 for m in zd.solution_metadata() if m.is_finite]

    result = SolveResult(solutions, run_id, where, int(zd.num_paths_hydrated()), zd)
    # top-level solves auto-declare their deliverable: "what were my solutions?"
    save("solutions [run %s]" % run_id, result,
         description="finite solutions, auto-declared by bertini.solve",
         directory=where)
    return result


def save(*args, description='', directory=None):
    """Save almost anything under a name: ``save(thing)`` or ``save(name, thing)``.

    A :class:`SolveResult` (or anything with ``.run_id`` and ``.solutions``) is declared
    as results with full provenance; any JSON-able value (dict, list, number, string) is
    recorded inline.  They land in the directory's ``results.json`` (pretty-printed,
    self-complete).  A nameless ``save(thing)`` is auto-named by timestamp;
    re-saving a name replaces it (newest wins).
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
    out.refresh_results()
    return name


def annotate(point, key, value, directory=None):
    """Attach metadata to a solution: ``annotate(sol, 'projection', 1.5)``.

    ``point`` is a :class:`Solution` (or anything with ``.provenance`` holding
    ``{'run', 'index'}``), or an explicit ``{'run': ..., 'index': ...}`` dict.  The
    annotation lands in the records beside the point it describes and renders into
    ``results.json``; re-annotating the same key replaces it (newest wins).  ``value``
    is any JSON-able thing.
    """
    provenance = getattr(point, 'provenance', None) or point
    if not isinstance(provenance, dict) or 'run' not in provenance or 'index' not in provenance:
        raise ValueError("annotate needs a point with provenance {'run', 'index'} "
                         "(a Solution from solve(), or an explicit dict)")
    out = _directory(directory)
    out.annotate(str(provenance['run']), int(provenance['index']), str(key),
                 _json.dumps(value))
    out.refresh_results()
    if hasattr(point, 'annotations'):
        point.annotations[str(key)] = value


def load(name=None, directory=None):
    """Load saved results by name -- the other half of :func:`save`.

    ``load("my favorites")`` returns that result (its points, annotations, provenance,
    or its inline value); ``load()`` returns the whole dict of everything saved, by
    name.  Reads the plain ``results.json``, so this works in any later session -- and
    the same file is readable without bertini at all.
    """
    path = _Path(directory if directory is not None else records_dir()) / 'results.json'
    if not path.exists():
        return {} if name is None else None
    document = _json.loads(path.read_text())
    # results.json is {"results": {...}, "runs": {...}}: the declared results plus
    # references to what constructed them.  load() serves the results section.
    everything = document.get('results', {})
    if name is None:
        return everything
    return everything.get(name)
