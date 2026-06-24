# This file is part of Bertini 2.
# 
# python/bertini/nag_algorithms/__init__.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/bertini/nag_algorithms/__init__.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/bertini/nag_algorithms/__init__.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#  silviana amethyst
#  UWEC
#  Spring 2023
# 





"""
nag_algorithms
"""


from bertini._pybertini import nag_algorithms as _pybnalag
from bertini._pybertini.nag_algorithms import *

# Ensure tracker classes have gained their .observers attribute (and PathDataCollector),
# which SolutionPathCollector relies on.  No cycle: bertini.tracking does not import this.
from .. import tracking as _tracking

# config structs gain update()/repr/to_dict/...; algorithm classes gain configure().
from ..config import _enhance_all, _enhance_owners
_enhance_all(_pybnalag)
_enhance_owners(_pybnalag)


# --- Slice: numpy-style row subsetting via [] and len(), on top of head()/tail()/rows() ---
#
# A Slice is a stack of linear forms; indexing/slicing it selects forms and returns a new (sub-)Slice
# over the same variables, so s[:k], s[-k:], s[[0, 2]], s[1] all compose.
def _slice_getitem(self, key):
    n = self.dimension()
    if isinstance(key, slice):
        return self.rows(list(range(*key.indices(n))))
    if isinstance(key, (list, tuple)):
        return self.rows([int(k) + n if int(k) < 0 else int(k) for k in key])
    k = int(key)
    return self.rows([k + n if k < 0 else k])


if hasattr(_pybnalag, 'Slice'):
    _pybnalag.Slice.__getitem__ = _slice_getitem
    _pybnalag.Slice.__len__ = lambda self: self.dimension()


# --- Bertini 1 / classic emission -----------------------------------------------------------
#
# A slice's linear forms are ordinary System functions (the classic writer expands the
# LinearFormsBlock), so a slice emits to a Bertini 1 input file via its standalone system.  A
# witness set emits its SQUARE system -- the witness set's system with the slice's forms appended
# (concatenate) -- which is what you would track in Bertini 1 for cross-validation.
def _slice_to_classic_input(self, **kwargs):
    """Emit this slice's linear forms as a Bertini 1 classic input file (see System.to_classic_input)."""
    return self.as_system().to_classic_input(**kwargs)


def _witness_system(self):
    """The square system whose isolated solutions are the witness points: this witness set's system
    with the slice's linear forms appended.  Requires the slice and system to share variables.

    Built by cloning the system (a copy sharing the variable nodes, so the witness set is left
    unmutated) and slice.add_to()-ing the clone.  add_to is homogenization-aware -- if the system was
    homogenized, it folds the slice's constant onto the homogenizing variable so the appended forms
    match the rest of the system -- which plain concatenate (a generic same-structure append) is not."""
    from bertini.system import clone
    sys = clone(self.get_system())
    self.get_slice().add_to(sys)
    return sys


def _witness_to_classic_input(self, **kwargs):
    """Emit the witness (square) system as a Bertini 1 classic input file (see System.to_classic_input)."""
    return self.witness_system().to_classic_input(**kwargs)


if hasattr(_pybnalag, 'Slice'):
    _pybnalag.Slice.to_classic_input = _slice_to_classic_input

for _ws_name in dir(_pybnalag):
    if _ws_name.startswith('WitnessSet'):
        _ws_cls = getattr(_pybnalag, _ws_name)
        if isinstance(_ws_cls, type):
            _ws_cls.witness_system = _witness_system
            _ws_cls.to_classic_input = _witness_to_classic_input


# --- to_dataframe: a ZeroDim solve as a pandas DataFrame -- the "database of solutions" ---
#
# One ROW per solution (per tracked path), columns = the coordinates (x0, x1, ...) followed by
# every field of the per-solution metadata.  Once it is a DataFrame, every category and its
# metadata is a one-line filter -- df[df.is_finite], df[df.is_real & ~df.is_singular],
# df[~df.is_finite] (at infinity) -- which is what the point/category accessors give you ad hoc.
#
# pandas is an OPTIONAL dependency (mirrors meta_observer's collector-as-DataFrame): the point
# accessors (all_solutions / finite_solutions / solution_metadata / ...) are the no-pandas path;
# to_dataframe raises a helpful ImportError if pandas is missing.

# The metadata columns, in a fixed order (so the DataFrame's shape is stable across solves and
# precision models).  Mirrors the fields exposed on SolutionMetaData (see zero_dim_export.hpp);
# the classification flags come first since they are what you filter on.
_SOLUTION_METADATA_FIELDS = (
    'path_index', 'solution_index',
    'is_finite', 'is_real', 'is_singular', 'multiplicity', 'multiplicity_representative',
    'condition_number', 'function_residual', 'newton_residual',
    'accuracy_estimate', 'accuracy_estimate_user_coords',
    'cycle_num', 'endgame_success', 'pre_endgame_success', 'final_time_used',
    'precision_changed', 'max_precision_used', 'time_of_first_prec_increase',
)


def _zerodim_to_dataframe(self, *, user_coords=True, omit_infinite=True, merge_multiplicities=True):
    """The solve as a pandas DataFrame -- one row per solution, the "database of solutions".

    Columns are ``solution`` -- the whole solution point, kept in a single cell -- then every
    per-solution metadata field (``is_finite``, ``is_real``, ``is_singular``, ``multiplicity``,
    ``condition_number``, ``endgame_success``, ``max_precision_used``, ...), and finally ``system``,
    a reference to the (target) system these solutions satisfy (so rows accumulated from several
    solves stay identifiable).  Each category is then a one-line filter, e.g.
    ``df[df.is_real & ~df.is_singular]`` (nonsingular real) or ``df[~df.is_finite]`` (at infinity).

    The ``solution`` cell is an independent copy of the solution vector (a numpy array of Python
    ``complex`` for a double solve, of :class:`bertini.multiprec.Complex` for a multiprecision one,
    so no precision is lost).  Coordinates are deliberately **not** exploded into ``x0, x1, ...``
    columns; split them yourself if you want them, e.g.
    ``df['x'] = [v[0] for v in df.solution]``.

    Parameters
    ----------
    user_coords : bool
        Coordinates in YOUR variables (dehomogenized; default), or the solver's internal
        homogenized on-patch coordinates when ``False`` -- the same choice as :meth:`all_solutions`.
    omit_infinite : bool
        Drop the endpoints not classified finite (``is_finite`` False) -- the at-infinity and
        failed paths.  ``True`` by default, so the frame holds just the genuine finite solutions;
        pass ``False`` to get every tracked path (their coordinate cells may be empty/NaN).
    merge_multiplicities : bool
        Collapse a multiplicity-``m`` solution -- which the solver returns as ``m`` coincident
        endpoints -- to its single representative row (``multiplicity`` still records ``m``).
        ``True`` by default.  Pass ``False`` to keep every endpoint, including the ``m-1``
        duplicate copies.  The grouping is the solver's own (the C++ clustering that computes
        multiplicity), read off ``multiplicity_representative``; this does not re-cluster.

    Returns
    -------
    pandas.DataFrame
        One row per solution; the ``solution`` cell is a copied vector whose elements are Python
        ``complex`` for a double-precision solve and :class:`bertini.multiprec.Complex` for a
        multiprecision one (kept native, so no precision is lost).

    Notes
    -----
    ``pandas`` is an optional dependency; this raises :class:`ImportError` if it is absent.  The
    point accessors -- :meth:`all_solutions`, :meth:`finite_solutions`, :meth:`solution_metadata`
    -- are the no-pandas path.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "ZeroDim.to_dataframe() needs the optional 'pandas' dependency (install it, e.g. "
            "`pip install pandas`).  Without pandas, use the point accessors instead: "
            "all_solutions(), finite_solutions(), solution_metadata()."
        ) from e

    points = self.all_solutions(user_coords)
    metadata = self.solution_metadata()
    system = self.target_system()       # the system these solutions satisfy (one shared reference)
    rows = []
    for i in range(min(len(points), len(metadata))):
        m = metadata[i]
        if omit_infinite and not m.is_finite:
            continue
        if merge_multiplicities and not m.multiplicity_representative:
            continue                    # a duplicate copy of an already-kept multiple solution
        # The whole solution lives in one cell.  .copy() snapshots the eigenpy vector's buffer --
        # essential, because indexing that vector lazily returns scalars that alias a reused
        # internal buffer, so storing the live vector (or its elements) and letting pandas read it
        # later would collapse every cell to one value.  The copy keeps the native element type, so
        # a multiprecision solve loses no precision.
        row = {'solution': points[i].copy()}
        for field in _SOLUTION_METADATA_FIELDS:
            row[field] = getattr(m, field)
        row['system'] = system          # a reference, so rows from different solves stay identifiable
        rows.append(row)

    return pd.DataFrame(rows)


def _attach_to_dataframe():
    """Give every bound ZeroDim class a to_dataframe method (idempotent)."""
    for name in dir(_pybnalag):
        if not name.startswith('ZeroDim'):
            continue
        cls = getattr(_pybnalag, name)
        if isinstance(cls, type) and not getattr(cls, '_b2_has_to_dataframe', False):
            cls.to_dataframe = _zerodim_to_dataframe
            cls._b2_has_to_dataframe = True


_attach_to_dataframe()


# --- ZeroDim: a friendly factory over the 18 bound ZeroDim<endgame x precision x start> classes ---

# Each bound class is named ZeroDim<Endgame><Precision>Precision<StartSystem>; rather than type
# ZeroDimCauchyFixedMultiplePrecisionTotalDegree, select the pieces by string.
_ZD_ENDGAMES = {
    'cauchy': 'Cauchy',
    'powerseries': 'PowerSeries', 'power_series': 'PowerSeries', 'pseg': 'PowerSeries',
}
_ZD_PRECISIONS = {
    'double': 'DoublePrecision', 'dbl': 'DoublePrecision', 'fixed_double': 'DoublePrecision',
    'multiple': 'FixedMultiplePrecision', 'mp': 'FixedMultiplePrecision',
    'fixed_multiple': 'FixedMultiplePrecision', 'multiprecision': 'FixedMultiplePrecision',
    'adaptive': 'AdaptivePrecision', 'amp': 'AdaptivePrecision',
}
_ZD_START_SYSTEMS = {
    'totaldegree': 'TotalDegree', 'total_degree': 'TotalDegree', 'td': 'TotalDegree',
    'mhom': 'MHomogeneous', 'mhomogeneous': 'MHomogeneous', 'multihomogeneous': 'MHomogeneous',
}


def _zd_select(value, table, kind):
    key = str(value).strip().lower().replace('-', '_')
    frag = table.get(key) or table.get(key.replace('_', ''))
    if frag is None:
        raise ValueError("ZeroDim: unknown {} {!r}; choose from {}"
                         .format(kind, value, sorted({k for k in table})))
    return frag


def _infer_start_system(system):
    """Pick the start system from the system's variable-group structure.

    Mirrors the C++ blackbox ``InferStartType`` (core/.../blackbox/switches_zerodim.hpp): a single
    affine variable group with no homogeneous/projective groups is the 1-homogeneous (total-degree)
    case; anything else -- two or more variable groups, or any homogeneous group -- is
    multihomogeneous, and total degree would be the wrong (over-counting) start system there.
    """
    if system.num_variable_groups() == 1 and system.num_hom_variable_groups() == 0:
        return 'totaldegree'
    return 'mhom'


def ZeroDim(system, *, endgame='cauchy', mptype='adaptive', startsystem='infer',
            precision=None):
    """Construct a zero-dim solver by name, with friendly defaults.

    ``ZeroDim(system)`` is the Cauchy endgame in adaptive precision with the start system **inferred
    from the system's variable-group structure** -- total degree for a single affine group,
    multihomogeneous otherwise -- so a multi-group (e.g. eigenvalue) system gets MHom automatically
    rather than an over-counting total-degree start.  Override any piece with a string::

        ZeroDim(system, endgame='cauchy', mptype='amp', startsystem='mhom')

    Parameters
    ----------
    system : the polynomial :class:`~bertini.System` to solve.
    endgame : ``'cauchy'`` (default) or ``'powerseries'``.
    mptype : the precision -- ``'double'``, ``'multiple'``, or ``'adaptive'`` (``'amp'``, the default).
    precision : an alias for ``mptype``; if given (not ``None``) it overrides ``mptype``.
    startsystem : ``'infer'`` (default -- choose from the variable-group structure, matching the
        C++ blackbox), or force it with ``'totaldegree'`` / ``'mhom'``.  To run from a homotopy you
        built yourself with given start points, use :func:`user_homotopy` / :func:`blend_homotopy`
        instead (their construction needs the homotopy and start points, not just a system).

    Returns a solver; call ``.solve()`` then ``.all_solutions()`` as for any zero-dim solver.

    Examples
    --------
    The default infers the start system; strings pick the rest::

        >>> import bertini
        >>> from bertini.nag_algorithm import ZeroDim
        >>> x = bertini.Variable('x')
        >>> sys = bertini.System()
        >>> sys.add_variable_group(bertini.VariableGroup([x]))
        >>> sys.add_function(x * x - 1)
        >>> type(ZeroDim(sys)).__name__
        'ZeroDimCauchyAdaptivePrecisionTotalDegree'
        >>> type(ZeroDim(sys, mptype='amp', startsystem='mhom')).__name__
        'ZeroDimCauchyAdaptivePrecisionMHomogeneous'
        >>> solver = ZeroDim(sys, mptype='adaptive')   # robust path
        >>> solver.solve()                             # doctest: +SKIP
        >>> solver.all_solutions()                         # doctest: +SKIP
    """
    if precision is not None:
        mptype = precision
    start_key = str(startsystem).strip().lower().replace('-', '_').replace('_', '')
    # user-homotopy can't be built from a system alone -- point at the right entry point.
    if start_key in ('user', 'userhomotopy'):
        raise ValueError(
            "ZeroDim does not build the user-homotopy solver (its construction needs a homotopy "
            "and start points, not just a system); build the homotopy with "
            "nag_algorithm.blend_homotopy / coefficient_parameter_homotopy and solve it with "
            "nag_algorithm.user_homotopy(homotopy, start_points, target).")
    if start_key == 'infer':
        startsystem = _infer_start_system(system)
    cls_name = ('ZeroDim'
                + _zd_select(endgame, _ZD_ENDGAMES, 'endgame')
                + _zd_select(mptype, _ZD_PRECISIONS, 'mptype')
                + _zd_select(startsystem, _ZD_START_SYSTEMS, 'startsystem'))
    return getattr(_pybnalag, cls_name)(system)


# --- user homotopy: run the zero-dim solver on a homotopy YOU built, from start points YOU have ---

_USER_HOMOTOPY_CLASSES = {
    ('double',   'cauchy'):      'ZeroDimCauchyDoublePrecisionUserHomotopy',
    ('double',   'powerseries'): 'ZeroDimPowerSeriesDoublePrecisionUserHomotopy',
    ('multiple', 'cauchy'):      'ZeroDimCauchyFixedMultiplePrecisionUserHomotopy',
    ('multiple', 'powerseries'): 'ZeroDimPowerSeriesFixedMultiplePrecisionUserHomotopy',
    ('adaptive', 'cauchy'):      'ZeroDimCauchyAdaptivePrecisionUserHomotopy',
    ('adaptive', 'powerseries'): 'ZeroDimPowerSeriesAdaptivePrecisionUserHomotopy',
}


class _UserHomotopySolver:
    """A thin holder around a user-homotopy ZeroDim solver.

    The underlying solver keeps *references* to the homotopy, the target system, and the start
    system, so this holder retains all three to keep them alive, and forwards every attribute
    and method (``solve``, ``solutions``, ``solution_metadata``, ``get_tracker``, ...) to the
    wrapped solver.
    """

    def __init__(self, solver, *kept_alive):
        object.__setattr__(self, '_solver', solver)
        object.__setattr__(self, '_kept_alive', kept_alive)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, '_solver'), name)


def user_homotopy(homotopy, start_points, target, *, precision='adaptive', endgame='cauchy'):
    """Run the zero-dim solver on a homotopy you constructed, from a list of start points you
    already have (e.g. the solutions of an earlier solve) -- the parameter-homotopy workflow.

    This reuses the entire zero-dim pipeline (pre-endgame tracking, the midpath check, the
    endgame, post-processing); it differs from the ``...TotalDegree`` / ``...MHomogeneous``
    entries only in that the homotopy and the start points are supplied, not generated.

    Parameters
    ----------
    homotopy : System
        The homotopy to track, with a path variable; tracked from the start time (default 1)
        down to 0.  Its t=1 slice must vanish at the given ``start_points``.
    start_points : iterable of vectors
        The start points (at the start time).  An earlier solve's ``solutions()`` works directly
        when the variable coordinates line up (e.g. an affine homotopy).
    target : System
        The system the solutions satisfy at t=0 -- used for dehomogenize / residual and for the
        solver's consistency check.  It must NOT have a path variable.
    precision : {'adaptive', 'double', 'multiple'}
        'adaptive' (default) is the robust path.
    endgame : {'cauchy', 'powerseries'}

    Returns a solver: call ``.solve()`` then ``.all_solutions()`` as for any zero-dim solver.
    """
    prec = {'double': 'double', 'multiple': 'multiple', 'fixed_multiple': 'multiple',
            'adaptive': 'adaptive'}.get(precision, precision)
    eg = {'cauchy': 'cauchy', 'powerseries': 'powerseries',
          'power_series': 'powerseries'}.get(endgame, endgame)
    try:
        cls_name = _USER_HOMOTOPY_CLASSES[(prec, eg)]
    except KeyError:
        raise ValueError(
            "user_homotopy: unknown (precision, endgame) = ({!r}, {!r}); "
            "precision in {{'adaptive','double','multiple'}}, endgame in {{'cauchy','powerseries'}}"
            .format(precision, endgame))
    solver_cls = getattr(_pybnalag, cls_name)
    # A frequent mix-up (issue #258): passing the start-point *solver* instead of its
    # start *points*.  A ZeroDim solver is not iterable, so list(start_points) below would
    # raise a cryptic "object is not iterable" naming an opaque class.  Catch it here and
    # say what to do.  (A list / numpy array / tuple of vectors has no .all_solutions/.get_tracker.)
    if hasattr(start_points, 'all_solutions') and hasattr(start_points, 'get_tracker'):
        raise TypeError(
            "user_homotopy: start_points must be the actual start *points* (an iterable of "
            "solution vectors), but a {} solver was passed.  Call its .solve() and then pass "
            "its .all_solutions():\n"
            "    start_solver.solve()\n"
            "    nag_algorithm.user_homotopy(homotopy, start_solver.all_solutions(), target)"
            .format(type(start_points).__name__))
    user_start = _pybnalag.UserStartSystem(target, list(start_points))
    solver = solver_cls(target, user_start, homotopy)
    return _UserHomotopySolver(solver, homotopy, target, user_start)


def coefficient_parameter_homotopy(target, generic, path_variable='t'):
    """Build a parameter homotopy interpolating two systems of the same shape.

    Returns H = (1 - t) * target + t * generic with ``t`` added as its path variable, so at t=1
    it is ``generic`` (whose solutions are your start points) and at t=0 it is ``target``.  Pair
    it with :func:`user_homotopy`: solve ``generic`` once, then reuse its solutions to move to
    ``target`` (and to any number of further targets that share ``generic``)::

        gen_solver = nag_algorithm.ZeroDim(generic, mptype='adaptive')
        gen_solver.solve()
        H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
        solver = nag_algorithm.user_homotopy(H, gen_solver.all_solutions(), target)
        solver.solve()

    ``target`` and ``generic`` must be built over the SAME variable objects (the interpolation
    combines their function trees).

    For robustness the generic system's coefficients should be *generic* (random complex), so the
    straight-line parameter path avoids the (measure-zero) singular locus.  This is the
    no-gamma-trick member of the family; if ``generic`` is *structured* (e.g. a products-of-linears
    start), it cannot be fused by System node arithmetic, so it is combined with a blend block --
    the same machinery as :func:`blend_homotopy`, but with the start coefficient fixed at 1.
    """
    # H = (1-t)*target + 1*t*generic.  Delegating to make_homotopy with gamma = the constant 1
    # (a) keeps the (1-t)/t semantics (no gamma trick) and (b) lets a structured-block ``generic``
    # be blended rather than SILENTLY DROPPED by System node arithmetic, which only combines the
    # polynomial block (see ADR-0020).
    from bertini.function_tree.symbol import Integer
    from bertini._pybertini import system as _system
    return _system.make_homotopy(target, generic, path_variable, Integer(1))


def blend_homotopy(target, start, *, path_variable='t', gamma=None):
    """Form the gamma-trick homotopy H = (1-t)*target + gamma*t*start for a start system you built.

    Unlike :func:`coefficient_parameter_homotopy` (node arithmetic, for two polynomial systems of the
    same shape), this also works when ``start`` carries a *structured evaluation block* -- e.g. a
    products-of-linears start system built with :func:`bertini.linalg.add_products_of_linears`.  Such
    a block cannot be fused by node arithmetic, so the two systems are combined with a blend block
    that evaluates whole Systems; this is the same construction the zero-dim solver uses internally
    for its generated (total-degree / multihomogeneous) start systems.

    Parameters
    ----------
    target : System
        The system whose solutions you want, reached at t=0.
    start : System
        A start system you authored, whose (known) solutions are the start points.  At t=1 the
        homotopy is ``gamma*start``, so those solutions are its roots.
    path_variable : str
        Name of the path variable t added to the homotopy (default ``'t'``).
    gamma : node or None
        The gamma coefficient.  ``None`` (default) draws a random rational gamma.  Pass an exact
        node (e.g. from :func:`bertini.linalg.coefficient`) off the real axis for a reproducible path.

    Returns
    -------
    System
        The homotopy; pair it with :func:`user_homotopy` and your start points to solve.
    """
    from bertini._pybertini import system as _system
    return _system.make_homotopy(target, start, path_variable, gamma)


def moving_homotopy(fixed, start_moving, end_moving, *, path_variable='t', gamma=None):
    """Form a homotopy that moves ONLY the moving rows, leaving the fixed system evaluated once.

    The regeneration / moving-slice homotopy:

        H = [ fixed's blocks ;  (1-t)*end_moving + gamma*t*start_moving ]

    ``fixed`` holds the equations that do not move -- the polynomial system and any *static* linear
    slices -- and stays as its own evaluation block(s); ``start_moving`` and ``end_moving`` hold just
    the rows that move (a linear slice that slides, or a products-of-linears that deforms into a
    polynomial), agreeing in function count and sharing ``fixed``'s variable structure.  Only the
    moving rows carry the path variable: the fixed blocks are evaluated once per point and contribute
    zero to ``dH/dt`` as the moving rows slide -- the fixed system is never re-evaluated or scaled by
    the path coefficient.

    At t=1 the moving rows are ``gamma*start_moving`` (so the start points are the roots of ``fixed``
    together with ``start_moving``); at t=0 they are ``end_moving``.  The fixed rows come first, then
    the moving rows; build the matching ``target`` for :func:`user_homotopy` as ``fixed`` concatenated
    with ``end_moving`` (e.g. via ``bertini.system.concatenate``).  ``gamma=None`` draws a random
    rational gamma.

    Returns the homotopy System; pair it with :func:`user_homotopy` and your start points to solve.
    """
    from bertini._pybertini import system as _system
    return _system.make_moving_homotopy(fixed, start_moving, end_moving, path_variable, gamma)


# --- SolutionPathCollector: collect every solution path of a whole solve, for plotting ---
#
# A two-level meta-observer.  Attach one to a ZeroDim solver; on each PathStarted it spins
# up a fresh tracking PathDataCollector, attaches it to the solver's tracker, and on the
# matching PathComplete harvests it into .series and detaches it.  Because the solver reuses
# one tracker for a path's main homotopy track AND its endgame sub-tracks, the per-path
# collector captures the *whole* journey to t -> 0 -- one clean series per solution path,
# endgame included (no start-time filtering needed).
class SolutionPathCollector(_pybnalag.observers.CustomObserver):
    """Collects each solution path of a ZeroDim solve into its own time series.

    Usage::

        a = SolutionPathCollector()
        solver.add_observer(a)
        solver.solve()
        for path in a.series:          # one tracking.PathDataCollector per solution path
            t, z = path.times(), path.points()
            ...
    """

    def __init__(self):
        super().__init__()
        self.series = []        # finished PathDataCollector per solution path, in completion order
        self._active = {}       # path_index -> (tracker, collector)

    def Observe(self, event):
        obs = _pybnalag.observers
        if isinstance(event, obs.PathStarted):
            tracker = event.solver().get_tracker()
            # tracker.observers is the precision-appropriate module (set in bertini.tracking)
            collector = tracker.observers.PathDataCollector()
            collector.path_index = event.path_index()
            tracker.add_observer(collector)
            self._active[event.path_index()] = (tracker, collector)
        elif isinstance(event, obs.PathComplete):
            entry = self._active.pop(event.path_index(), None)
            if entry is not None:
                tracker, collector = entry
                tracker.remove_observer(collector)
                self.series.append(collector)


_pybnalag.observers.SolutionPathCollector = SolutionPathCollector


__all__ = dir(_pybnalag)
__all__.append('ZeroDim')
__all__.append('user_homotopy')
__all__.append('coefficient_parameter_homotopy')
__all__.append('moving_homotopy')
__all__.append('blend_homotopy')
__all__.append('SolutionPathCollector')


# DoublePrecisionTotalDegree = bertini._pybertini.nag_algorithms.ZeroDimCauchyDoublePrecisionTotalDegree