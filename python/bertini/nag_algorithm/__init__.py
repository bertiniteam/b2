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


# --- Slice: a Python sequence of linear forms ----------------------------------------------
#
# A Slice is a stack of linear forms; treat it as a Python sequence, following list semantics:
#   * an integer index selects an ELEMENT -- the i-th form's coefficient VECTOR (a 1-D array);
#     "a single linear form is a vector".  Iterating the slice yields these form vectors.
#   * a Python slice (s[i:j]) selects a SUB-COLLECTION -- a new (sub-)Slice (same type).  A list /
#     tuple of indices likewise selects a sub-Slice (a numpy-flavoured multi-select).
# The whole coefficient matrix is `coefficients()` (always 2-D); see the coefficients wrapper below.
#
# The vector-vs-matrix distinction is thus explicit (element index vs collection slice / matrix
# accessor), never an emergent property of the form count -- which is exactly the eigenpy 1-D
# "collapse" footgun we are insulating callers from.  See docs/adr for the eigenpy shape rule.

# eigenpy collapses a 1-row Eigen matrix to a 1-D numpy array, so the bound coefficients() is 1-D for
# a single-form slice and 2-D otherwise -- a data-dependent shape.  Wrap it to ALWAYS return 2-D,
# reshaped to (num_forms, num_variables+1) using the dimensions we know in C++ (orientation-correct).
_slice_coefficients_native = _pybnalag.Slice.coefficients if hasattr(_pybnalag, 'Slice') else None


def _slice_coefficients(self):
    """The augmented coefficient matrix of the slice's linear forms -- **always 2-D**.

    Shape ``(num_forms, num_variables + 1)``: one row per linear form, the trailing column carrying
    each form's constant term.  This holds even for a single-form slice, where the underlying
    binding library would otherwise hand back a 1-D array (eigenpy collapses a one-row matrix).  The
    2-D shape is part of *this* accessor's contract -- ours, not the binding's -- so it is stable
    across binding libraries; see ``docs/adr/0033``.  For one form's coefficient *vector*, index an
    element: ``slice[i]``.

    Examples
    --------
    >>> import bertini                                                   # doctest: +SKIP
    >>> import bertini                                                   # doctest: +SKIP
    >>> x, y = bertini.Variable('x'), bertini.Variable('y')             # doctest: +SKIP
    >>> bertini.Slice.from_coefficients([[2, 3, 1]], [x, y]).coefficients().shape  # doctest: +SKIP
    (1, 3)
    """
    import numpy as np
    raw = np.asarray(_slice_coefficients_native(self))
    return raw.reshape(self.dimension(), self.num_variables() + 1)


def _slice_getitem(self, key):
    """Index a slice like a Python sequence of its linear forms.

    * ``slice[i]`` (an integer) -> the i-th form's coefficient **vector**, a 1-D array of length
      ``num_variables + 1`` (trailing entry = constant term).  Iterating (``for form in slice``)
      yields these vectors.  "A single linear form is a vector."
    * ``slice[i:j]`` (a Python slice) -> a **sub-Slice**, a new Slice of those forms.  A list/tuple
      of indices (``slice[[0, 2]]``) likewise -> a sub-Slice.

    So the vector view and the matrix view are *named* (element index vs slice / :meth:`coefficients`),
    never inferred from the form count.  See ``docs/adr/0033``.

    Examples
    --------
    >>> import bertini                                                   # doctest: +SKIP
    >>> import bertini                                                   # doctest: +SKIP
    >>> x, y = bertini.Variable('x'), bertini.Variable('y')             # doctest: +SKIP
    >>> s = bertini.Slice.from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])  # doctest: +SKIP
    >>> s[0]                       # an ELEMENT: the first form's vector  # doctest: +SKIP
    array([2, 3, 1], dtype=object)
    >>> s[:1].dimension()          # a SUB-COLLECTION: a one-form Slice   # doctest: +SKIP
    1
    """
    n = self.dimension()
    if isinstance(key, slice):
        return self.rows(list(range(*key.indices(n))))                 # sub-Slice (a sub-collection)
    if isinstance(key, (list, tuple)):
        return self.rows([int(k) + n if int(k) < 0 else int(k) for k in key])   # sub-Slice
    k = int(key)
    if k < 0:
        k += n
    if not (0 <= k < n):
        raise IndexError("slice form index {!r} out of range [0, {})".format(key, n))
    return self.coefficients()[k]                                      # element: the form's vector


_SLICE_DOC = """A linear slice: a stack of linear forms M [x ; 1] that cuts a positive-dimensional
component down to witness points.  The linear part of a witness set.

A slice is a Python **sequence of its linear forms**:

* ``len(slice)`` -- the number of forms (the slice's dimension).
* ``slice[i]`` (integer) -- the i-th form's coefficient **vector** (1-D, length num_variables+1).
  Iterating yields these vectors.  A single linear form is a vector.
* ``slice[i:j]`` / ``slice[[i, j]]`` -- a **sub-Slice** (a sub-collection of forms).
* ``slice.coefficients()`` -- the whole augmented coefficient matrix, **always 2-D**
  ``(num_forms, num_variables+1)``.

The vector view vs the matrix view is named (element index vs slice / coefficients), never inferred
from the form count -- so it is stable regardless of the binding library's shape conventions
(docs/adr/0033).

A slice does NOT own homogenization -- the system does.  ``slice.add_to(system)`` appends the forms
to a system (folding the constant onto the homogenizing variable if the system was homogenized);
``slice.as_system()`` returns a standalone System of just the forms.  Build slices with
``Slice.random_complex`` / ``Slice.random_real`` / ``Slice.from_coefficients`` (or
``bertini.bertini.Slice.from_coefficients`` for exact numpy/list coefficients).
"""

if hasattr(_pybnalag, 'Slice'):
    _pybnalag.Slice.coefficients = _slice_coefficients
    _pybnalag.Slice.__getitem__ = _slice_getitem
    _pybnalag.Slice.__len__ = lambda self: self.dimension()
    try:
        _pybnalag.Slice.__doc__ = _SLICE_DOC
    except (AttributeError, TypeError):
        pass   # some Boost.Python builds make the class docstring read-only; the methods carry docs


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


def _plural(n, noun):
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def _witness_repr(self):
    """A concise summary: the component's dimension and degree, and a one-line description of the
    system and slice it carries.  (The slice and system have their own detailed reprs.)"""
    try:
        consistent = 'consistent' if self.is_consistent() else 'INCONSISTENT'
    except Exception:
        consistent = 'consistency unknown'

    lines = ["WitnessSet -- dimension {}, degree {} ({})".format(
        self.dimension(), self.degree(), consistent)]

    try:
        s = self.get_slice()
        lines.append("  slice:  {} on {}".format(
            _plural(s.dimension(), 'linear form'), _plural(s.num_variables(), 'variable')))
    except Exception:
        lines.append("  slice:  <none>")

    try:
        sysm = self.get_system()
        lines.append("  system: {} on {}".format(
            _plural(sysm.num_functions(), 'function'), _plural(sysm.num_variables(), 'variable')))
    except Exception:
        lines.append("  system: <none>")

    return "\n".join(lines)


for _ws_name in dir(_pybnalag):
    if _ws_name.startswith('WitnessSet'):
        _ws_cls = getattr(_pybnalag, _ws_name)
        if isinstance(_ws_cls, type):
            _ws_cls.witness_system = _witness_system
            _ws_cls.to_classic_input = _witness_to_classic_input
            _ws_cls.__repr__ = _witness_repr
            _ws_cls.__str__ = _witness_repr


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
    'is_finite', 'is_real', 'is_singular', 'is_nonsolution', 'crossing_unresolved',
    'multiplicity', 'multiplicity_representative',
    'condition_number', 'function_residual', 'newton_residual',
    'accuracy_estimate', 'accuracy_estimate_user_coords',
    'precision_digits', 'accuracy_digits',
    'cycle_num', 'endgame_success_code', 'pre_endgame_success_code', 'final_time_used',
    'precision_changed', 'max_precision_used', 'time_of_first_prec_increase',
    'path_time_seconds',
    'num_successful_steps', 'num_failed_steps', 'latest_path_point', 'wall_clock_limit_seconds',
)


def _zerodim_to_dataframe(self, *, user_coords=True, omit_infinite=True, merge_multiplicities=True,
                          group=None):
    """The solve as a pandas DataFrame -- one row per solution, the "database of solutions".

    Columns are ``solution`` -- the whole solution point, kept in a single cell -- then every
    per-solution metadata field (``is_finite``, ``is_real``, ``is_singular``, ``multiplicity``,
    ``condition_number``, ``endgame_success_code``, ``max_precision_used``, ...), and finally ``system``,
    a reference to the (target) system these solutions satisfy (so rows accumulated from several
    solves stay identifiable).  Each category is then a one-line filter, e.g.
    ``df[df.is_real & ~df.is_singular]`` (nonsingular real) or ``df[~df.is_finite]`` (at infinity).

    The ``solution`` cell is an independent copy of the solution vector (a numpy array of Python
    ``complex`` for a double solve, of :class:`bertini.complex_mp` for a multiprecision one,
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
    group : VariableGroup or int, optional
        Project each solution onto one variable group -- the ``solution`` cell then holds only that
        group's coordinates (Cluster G).  Pass the ``VariableGroup`` object or its 0-based FIFO index.
        ``None`` (default) keeps the whole point.

    Returns
    -------
    pandas.DataFrame
        One row per solution; the ``solution`` cell is a copied vector whose elements are Python
        ``complex`` for a double-precision solve and :class:`bertini.complex_mp` for a
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
        sol = points[i].copy()
        if group is not None:
            sol = system.coordinates_of(sol, group)   # project onto one variable group (Cluster G)
        row = {'solution': sol}
        for field in _SOLUTION_METADATA_FIELDS:
            row[field] = getattr(m, field)
        row['system'] = system          # a reference, so rows from different solves stay identifiable
        rows.append(row)

    return pd.DataFrame(rows)


def _attach_to_dataframe():
    """Give every bound solver class (ZeroDimSolver*, HomotopySolver*) a to_dataframe method (idempotent)."""
    for name in dir(_pybnalag):
        if not name.startswith(('ZeroDimSolver', 'HomotopySolver')):
            continue
        cls = getattr(_pybnalag, name)
        if isinstance(cls, type) and not getattr(cls, '_b2_has_to_dataframe', False):
            cls.to_dataframe = _zerodim_to_dataframe
            cls._b2_has_to_dataframe = True


_attach_to_dataframe()


# --- group= : project each returned solution onto one variable group (Cluster G) ------------
#
# The reusable primitive is C++ (System.coordinates_of); here we just surface a `group=` keyword on
# every point getter so `solver.solutions(group=xyz)` returns only that group's coordinates.  The
# projection goes through the solver's target_system (whose FIFO layout the user-coord points follow).
_POINT_GETTERS_FOR_GROUP = (
    'all_solutions', 'solutions', 'finite_solutions', 'real_solutions',
    'nonsingular_solutions', 'singular_solutions', 'infinite_solutions', 'nonsolutions',
)


def _make_group_getter(native):
    def getter(self, *args, group=None, **kwargs):
        points = native(self, *args, **kwargs)
        if group is None:
            return points
        sysm = self.target_system()
        return [sysm.coordinates_of(p, group) for p in points]
    getter.__name__ = getattr(native, '__name__', 'solutions')
    getter.__doc__ = (getattr(native, '__doc__', '') or '') + (
        "\n\ngroup=<VariableGroup or 0-based FIFO index>: project each returned point onto that "
        "variable group -- return only its coordinates (Cluster G).  Default None keeps the whole point.")
    return getter


def _attach_group_kwarg():
    """Give every solver class a group= keyword on its point getters (idempotent)."""
    for name in dir(_pybnalag):
        if not name.startswith(('ZeroDimSolver', 'HomotopySolver')):
            continue
        cls = getattr(_pybnalag, name)
        if not isinstance(cls, type) or getattr(cls, '_b2_has_group_kwarg', False):
            continue
        for getter_name in _POINT_GETTERS_FOR_GROUP:
            native = getattr(cls, getter_name, None)
            if native is not None:
                setattr(cls, getter_name, _make_group_getter(native))
        cls._b2_has_group_kwarg = True


_attach_group_kwarg()


# --- solve() returns a SolveResult; result() re-derives it -----------------------------------
#
# The solver records ITSELF (record_to / the ambient directory), so it can hand back the same
# records-aware SolveResult that bertini.solve returns -- no dependence on the records-layer
# orchestration.  This makes solver.solve() and bertini.solve() return the same type (removing the
# "the bare solver returns nothing" inconsistency), and result() re-derives it if the caller
# dropped the return value.
def _solver_result(self):
    """This solve's :class:`~bertini.records.SolveResult`: the finite solutions plus the records
    ticket (run id, directory, recall count), read from the solver's own recorded state.

    A bare ``solver.solve()`` already returns this; ``result()`` re-derives it (call it after
    ``solve()``) if you did not keep the return value.
    """
    from bertini.records import SolveResult
    return SolveResult.from_solver(self)


def _describe_solver_class(cls_name):
    """A friendly label from a bound solver class name, e.g.
    'ZeroDimSolverCauchyAdaptivePrecision' -> 'ZeroDimSolver[cauchy, adaptive precision]'."""
    for kind in ('ZeroDimSolver', 'HomotopySolver'):
        if cls_name.startswith(kind):
            rest = cls_name[len(kind):]
            eg = ('cauchy' if 'Cauchy' in rest else
                  'power-series' if 'PowerSeries' in rest else '?')
            prec = ('double' if 'Double' in rest else
                    'fixed-multiple' if 'FixedMultiple' in rest else
                    'adaptive' if 'Adaptive' in rest else '?')
            return '{}[{}, {} precision]'.format(kind, eg, prec)
    return cls_name


def _solver_repr(self):
    """A readable one-line summary: the solver kind, and once solved its solution tally.

    Replaces the useless default ``<...object at 0x...>``.  Before solving: the kind plus a nudge
    to call ``.solve()``.  After: finite/real/singular counts out of the paths tracked, and the
    records run id when the solve recorded.
    """
    label = _describe_solver_class(type(self).__name__)
    try:
        md = self.solution_metadata()
    except Exception:
        return '<{}>'.format(label)
    if not len(md):
        return '<{}: not yet solved -- call .solve()>'.format(label)
    # count DISTINCT finite solutions (multiplicity representatives), matching the default
    # merge_multiplicities view of finite_solutions(); len(md) is the raw path count.
    reps = [m for m in md if m.is_finite and m.multiplicity_representative]
    n_finite = len(reps)
    n_real = sum(1 for m in reps if m.is_real)
    n_sing = sum(1 for m in reps if m.is_singular)
    run = self.records_run_id()
    tail = ', run {}'.format(run) if run else ''
    return ('<{}: {} finite solution{} ({} real, {} singular) of {} path{}{}>'
            .format(label, n_finite, '' if n_finite == 1 else 's',
                    n_real, n_sing, len(md), '' if len(md) == 1 else 's', tail))


def _make_solve_returning_result(native_solve):
    def solve(self, communicator=None, settings=None, **field_settings):
        # settings= : a dict of config fields; **field_settings: the same by keyword.  Both are
        # applied via set() before solving -- collapsing make/set/solve to one line.  communicator
        # and settings are reserved names (MPI is a solve-time opt-in), never config fields.
        combined = dict(settings or {}, **field_settings)
        if combined:
            self.set(**combined)
        native_solve(self, communicator)
        return self.result()
    solve.__name__ = 'solve'
    solve.__doc__ = ((getattr(native_solve, '__doc__', '') or '') +
        "\n\nSettings: pass config fields either as a dict (settings={'final_tolerance': 1e-13}) or by "
        "keyword (solve(final_tolerance=1e-13)); each is applied via set() before solving, so "
        "make/set/solve becomes a single call.  communicator= is reserved for MPI (never a setting).\n\n"
        "Returns a bertini.records.SolveResult -- the finite solutions plus this solve's records "
        "ticket (run id, directory, recall count), read from the solver's own records (the solver "
        "records itself, on by default).  Drop it freely: solver.result() re-derives it, and the "
        "records hold the truth.")
    return solve


def _attach_result_and_solve():
    """Make every solver's solve() return a SolveResult and give it result() (idempotent)."""
    for name in dir(_pybnalag):
        if not name.startswith(('ZeroDimSolver', 'HomotopySolver')):
            continue
        cls = getattr(_pybnalag, name)
        if not isinstance(cls, type) or getattr(cls, '_b2_has_result', False):
            continue
        cls.result = _solver_result
        cls.solve = _make_solve_returning_result(cls.solve)
        cls.__repr__ = _solver_repr
        cls._b2_has_result = True


_attach_result_and_solve()


# --- ZeroDimSolver: a friendly factory over the bound ZeroDimSolver<endgame x precision> classes ---

# Each bound solver class is named ZeroDimSolver<Endgame><Precision> (the start system is NO
# longer part of the type -- it is chosen at construction).  Select endgame + precision by string to
# pick the class, then pick the start system separately (default total_degree_binomial, else a factory).
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
# string -> the StartSystemType enum value name on _pybnalag (the C++ blackbox::type::Start).
_ZD_START_SYSTEMS = {
    'totaldegreebinomial': 'total_degree_binomial', 'binomial': 'total_degree_binomial', 'tdb': 'total_degree_binomial',
    'totaldegreelinearproduct': 'total_degree_linear_product', 'linearproduct': 'total_degree_linear_product', 'tdlp': 'total_degree_linear_product',
    'mhom': 'mhomogeneous', 'mhomogeneous': 'mhomogeneous', 'multihomogeneous': 'mhomogeneous',
}


def _zd_select(value, table, kind):
    key = str(value).strip().lower().replace('-', '_')
    frag = table.get(key) or table.get(key.replace('_', ''))
    if frag is None:
        raise ValueError("ZeroDimSolver: unknown {} {!r}; choose from {}"
                         .format(kind, value, sorted({k for k in table})))
    return frag


def _infer_start_system(system):
    """Pick the start system from the system's variable-group structure.

    Mirrors the C++ blackbox ``InferStartType`` (core/.../blackbox/switches_zerodim.hpp): a single
    affine variable group with no homogeneous/projective groups is the 1-homogeneous case; anything
    else -- two or more variable groups, or any homogeneous group -- is multihomogeneous (a
    total-degree start would be the wrong, over-counting choice there).

    The single-affine-group case infers ``total_degree_binomial`` -- the default: it is as reliable
    as the linear-product total-degree start and faster, so it stays the default for 1-homogeneous /
    1-affine systems.  See core/.../nag_algorithms/common/policies.hpp.
    """
    if system.num_variable_groups() == 1 and system.num_hom_variable_groups() == 0:
        return 'binomial'
    return 'mhom'


def _precision_model(mptype, precision):
    """Resolve the precision *model* (mptype) and apply an integer *digit count* (precision).

    ``precision`` is the integer number of digits; it is applied via ``bertini.default_precision``
    at construction (the documented "set precision, then build" step, which the fixed-multiple and
    adaptive trackers read at construction).  ``mptype`` is the precision MODEL
    (``'double'`` / ``'multiple'`` / ``'adaptive'``).

    A **string** ``precision`` is the old, confusing model-selector alias -- honored as ``mptype``
    with a ``DeprecationWarning`` for one release.  Returns the mptype string to build with.
    """
    import warnings
    if isinstance(precision, str):
        warnings.warn(
            "precision={0!r} as a precision-MODEL is deprecated and will be removed: pass "
            "mptype={0!r} for the model, and reserve precision= for an integer number of digits."
            .format(precision), DeprecationWarning, stacklevel=3)
        return precision
    if precision is not None:
        import bertini
        bertini.default_precision(int(precision))     # the "set precision, then build" step
    return mptype


def ZeroDimSolver(system, *, endgame='powerseries', mptype='adaptive', startsystem='infer',
                  precision=None, settings=None, **field_settings):
    """Construct a zero-dim solver by name, with friendly defaults.

    ``ZeroDimSolver(system)`` is the Cauchy endgame in adaptive precision with the start system
    **inferred from the system's variable-group structure** -- total degree for a single affine group,
    multihomogeneous otherwise -- so a multi-group (e.g. eigenvalue) system gets MHom automatically
    rather than an over-counting total-degree start.  Override any piece with a string::

        ZeroDimSolver(system, endgame='cauchy', mptype='amp', startsystem='mhom')

    Config settings may be applied at construction (via ``set``) so no separate ``set()`` is needed
    before solving -- either as a dict (``ZeroDimSolver(sys, settings={'final_tolerance': 1e-13})``)
    or by keyword (``ZeroDimSolver(sys, final_tolerance=1e-13)``); the same go to ``solve(...)`` too.

    Parameters
    ----------
    system : the polynomial :class:`~bertini.System` to solve.
    endgame : ``'powerseries'`` (default, matching Bertini 1; ADR-0058) or ``'cauchy'``.
    mptype : the precision MODEL -- ``'double'``, ``'multiple'``, or ``'adaptive'`` (``'amp'``, the default).
    precision : the number of DIGITS (an ``int``), applied via ``bertini.default_precision`` at
        construction -- meaningful for ``'multiple'``/``'adaptive'`` (``'double'`` is always 16).  A
        *string* here is the deprecated old spelling of ``mptype`` and warns.
    startsystem : ``'infer'`` (default -- choose from the variable-group structure, matching the
        C++ blackbox), or force it with ``'binomial'`` / ``'linearproduct'`` / ``'mhom'``.  To run from a homotopy you
        built yourself with given start points, use :class:`HomotopySolver` / :func:`straight_line_homotopy`
        instead (their construction needs the homotopy and start points, not just a system).

    Returns a solver; call ``.solve()`` then ``.all_solutions()`` as for any zero-dim solver.

    Examples
    --------
    The default infers the start system; strings pick the rest::

        >>> import bertini
        >>> from bertini.nag_algorithm import ZeroDimSolver
        >>> x = bertini.Variable('x')
        >>> sys = bertini.System()
        >>> sys.add_variable_group(bertini.VariableGroup([x]))
        >>> sys.add_function(x * x - 1)
        >>> type(ZeroDimSolver(sys)).__name__
        'ZeroDimSolverCauchyAdaptivePrecision'
        >>> type(ZeroDimSolver(sys, mptype='amp', startsystem='mhom')).__name__
        'ZeroDimSolverCauchyAdaptivePrecision'
        >>> solver = ZeroDimSolver(sys, mptype='adaptive')   # robust path
        >>> solver.solve()                             # doctest: +SKIP
        >>> solver.all_solutions()                         # doctest: +SKIP
    """
    mptype = _precision_model(mptype, precision)
    start_key = str(startsystem).strip().lower().replace('-', '_').replace('_', '')
    # user-homotopy can't be built from a system alone -- point at the right entry point.
    if start_key in ('user', 'userhomotopy'):
        raise ValueError(
            "ZeroDimSolver does not build the homotopy solver (its construction needs a homotopy "
            "and start points, not just a system); build the homotopy with "
            "nag_algorithm.straight_line_homotopy and solve it with "
            "nag_algorithm.HomotopySolver(homotopy, start_points, target).")
    if start_key == 'infer':
        startsystem = _infer_start_system(system)
    # the solver class encodes only endgame + precision now; the start system is a constructor choice.
    cls_name = ('ZeroDimSolver'
                + _zd_select(endgame, _ZD_ENDGAMES, 'endgame')
                + _zd_select(mptype, _ZD_PRECISIONS, 'mptype'))
    cls = getattr(_pybnalag, cls_name)
    start_enum = getattr(_pybnalag.StartSystemType,
                         _zd_select(startsystem, _ZD_START_SYSTEMS, 'startsystem'))
    # Every start system is built through its explicit factory.  (Previously total_degree
    # short-circuited to the default constructor, which silently substituted the *default* start
    # system -- the facade quirk that made startsystem='totaldegree' secretly solve the default.)
    solver = cls(system, _pybnalag.start_system_factory(start_enum))
    # config fields as a dict (settings={...}) and/or by keyword, applied at construction so a one-line
    # ZeroDimSolver(sys, final_tolerance=1e-13) needs no separate set() before solve().
    combined = dict(settings or {}, **field_settings)
    if combined:
        solver.set(**combined)
    return solver


# --- HomotopySolver: track a homotopy YOU built, from start points YOU have ----------------------
#
# The continuation primitive: given a homotopy with a path variable and a list of start points,
# track each path through the endgame and classify the endpoints.  This is the SAME machinery the
# zero-dim solve uses, with the homotopy + start points supplied rather than generated.

_HOMOTOPY_SOLVER_CLASSES = {
    ('double',   'cauchy'):      'HomotopySolverCauchyDoublePrecision',
    ('double',   'powerseries'): 'HomotopySolverPowerSeriesDoublePrecision',
    ('multiple', 'cauchy'):      'HomotopySolverCauchyFixedMultiplePrecision',
    ('multiple', 'powerseries'): 'HomotopySolverPowerSeriesFixedMultiplePrecision',
    ('adaptive', 'cauchy'):      'HomotopySolverCauchyAdaptivePrecision',
    ('adaptive', 'powerseries'): 'HomotopySolverPowerSeriesAdaptivePrecision',
}


class _HomotopySolverHolder:
    """A thin holder around a bound HomotopySolver.

    The underlying solver keeps *references* to the homotopy, the target system, and the start
    system, so this holder retains all three to keep them alive, and forwards every attribute
    and method (``solve``, ``all_solutions``, ``solution_metadata``, ``get_tracker``, ...) to the
    wrapped solver.
    """

    def __init__(self, solver, *kept_alive):
        object.__setattr__(self, '_solver', solver)
        object.__setattr__(self, '_kept_alive', kept_alive)

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, '_solver'), name)

    def __repr__(self):
        # special methods bypass __getattr__, so delegate the friendly repr explicitly
        return repr(object.__getattribute__(self, '_solver'))


def _coerce_start_point_coordinate(entry, which_point, which_coordinate):
    """One start-point coordinate, as a ``complex_mp`` the native layer accepts.

    Start points are *transported values* (the tracker refines them), so anything that IS a
    number is welcome: multiprecision scalars, Python/numpy numerics, and CONSTANT symbolic
    nodes (evaluated).  An expression still containing variables is not a number, and gets a
    math error saying so; anything else gets a TypeError naming the accepted kinds.
    """
    from bertini._pybertini.multiprec import complex_mp as _complex_mp, real_mp as _real_mp
    from bertini._pybertini.function_tree import AbstractNode as _AbstractNode

    if isinstance(entry, _complex_mp):
        return entry
    if isinstance(entry, _real_mp):
        return _complex_mp(entry)
    if isinstance(entry, _AbstractNode):
        vars_inside = [str(v) for v in entry.variables()]
        if vars_inside:
            raise ValueError(
                "start point {} coordinate {}: a start point must be a NUMBER, but this "
                "coordinate is a symbolic expression still containing the variable(s) {} -- "
                "evaluate or substitute it down to a constant first"
                .format(which_point, which_coordinate, vars_inside))
        return entry.eval()
    if isinstance(entry, bool):
        # bool is an int subclass; a True/False coordinate is a bug, not a number.
        raise TypeError(
            "start point {} coordinate {}: got a bool; start-point coordinates must be numbers"
            .format(which_point, which_coordinate))
    if isinstance(entry, int):
        return _complex_mp(entry)
    try:
        z = complex(entry)  # Python float/complex and every numpy scalar kind land here
    except (TypeError, ValueError):
        raise TypeError(
            "start point {} coordinate {}: cannot interpret a {} as a number.  Start-point "
            "coordinates may be bertini multiprecision scalars (complex_mp / real_mp), Python "
            "or numpy numbers, or constant symbolic nodes (e.g. bertini.symbolics.Complex)"
            .format(which_point, which_coordinate, type(entry).__name__)) from None
    return _complex_mp(z.real, z.imag)


def _coerce_start_points(start_points):
    """The user's start points, normalized to vectors of ``complex_mp``.

    This is the tolerant seam between "whatever numbers the user has" and the native solver,
    which needs multiprecision vectors: each point may be a list/tuple/numpy array (any
    numeric dtype, including object arrays of mp values or constant symbolic nodes), and each
    coordinate goes through :func:`_coerce_start_point_coordinate`.  Issue #347.
    """
    import numpy as np
    pts = []
    for k, pt in enumerate(start_points):
        if isinstance(pt, (str, bytes)) or not hasattr(pt, '__iter__'):
            raise TypeError(
                "start point {}: each start point must be a vector of coordinates (a list, "
                "tuple, or numpy array), got a {}".format(k, type(pt).__name__))
        coerced = [_coerce_start_point_coordinate(e, k, j) for j, e in enumerate(pt)]
        # no dtype= here: complex_mp registers its own numpy dtype, and inference finds it;
        # an explicit dtype=object array is NOT accepted by the native converter.
        pts.append(np.array(coerced))
    return pts


def HomotopySolver(homotopy, start_points, target, *, mptype='adaptive', precision=None, endgame='powerseries', amp_config=None):
    """Track a homotopy you constructed, from a list of start points you already have (e.g. the
    solutions of an earlier solve) -- the continuation primitive (parameter-homotopy workflow).

    This reuses the entire tracking pipeline (pre-endgame tracking, the midpath check, the
    endgame, post-processing); it differs from :func:`ZeroDimSolver` only in that the homotopy and
    the start points are supplied, not generated.

    Parameters
    ----------
    homotopy : System
        The homotopy to track, with a path variable; tracked from the start time (default 1)
        down to 0.  Its t=1 slice must vanish at the given ``start_points``.
    start_points : iterable of vectors
        The start points (at the start time).  An earlier solve's ``all_solutions()`` works directly
        when the variable coordinates line up (e.g. an affine homotopy).  Each point may be a list,
        tuple, or numpy array; coordinates may be multiprecision scalars (``complex_mp`` /
        ``real_mp``), Python or numpy numbers, or constant symbolic nodes -- they are all coerced
        to multiprecision for you.  (Start points are transported values, refined by the tracker,
        so double-precision input is fine here.)
    target : System
        The system the solutions satisfy at t=0 -- used for dehomogenize / residual and for the
        solver's consistency check.  It must NOT have a path variable.
    mptype : {'adaptive', 'double', 'multiple'}
        The precision MODEL; 'adaptive' (default) is the robust path.  Adaptive precision needs a
        degree bound, so it is refused for a homotopy that is not polynomial unless you supply
        ``amp_config``; such a homotopy tracks at fixed precision without any of this.
    amp_config : AMPConfig, optional
        An adaptive-precision configuration to install on the tracker instead of the one derived
        from the system.  This is how an analytic homotopy gets adaptive precision: the derivation
        needs a degree the system does not have, so you choose ``jacobian_eval_error_bound`` and
        ``function_eval_error_bound`` yourself.  Ignored unless ``mptype`` is 'adaptive'.
    precision : int, optional
        The number of DIGITS, applied via ``bertini.default_precision`` at construction.  A string
        here is the deprecated old spelling of ``mptype`` and warns.
    endgame : {'cauchy', 'powerseries'}

    Returns a solver: call ``.solve()`` then ``.all_solutions()`` as for any solver.

    Whether the homotopy is projective is settled where it was *built* -- see
    :func:`straight_line_homotopy` and :func:`moving_homotopy`, which projectivize by default.
    This function only tracks what it is given.  Start points written in your own (affine)
    coordinates are lifted onto the homotopy's patch for you, so an earlier affine solve's
    solutions feed a projective homotopy directly.
    """
    mptype = _precision_model(mptype, precision)
    prec = {'double': 'double', 'multiple': 'multiple', 'fixed_multiple': 'multiple',
            'adaptive': 'adaptive'}.get(mptype, mptype)
    eg = {'cauchy': 'cauchy', 'powerseries': 'powerseries',
          'power_series': 'powerseries'}.get(endgame, endgame)
    try:
        cls_name = _HOMOTOPY_SOLVER_CLASSES[(prec, eg)]
    except KeyError:
        raise ValueError(
            "HomotopySolver: unknown (mptype, endgame) = ({!r}, {!r}); "
            "mptype in {{'adaptive','double','multiple'}}, endgame in {{'cauchy','powerseries'}}"
            .format(mptype, endgame))
    solver_cls = getattr(_pybnalag, cls_name)
    # A frequent mix-up (issue #258): passing the start-point *solver* instead of its
    # start *points*.  A solver is not iterable, so list(start_points) below would
    # raise a cryptic "object is not iterable" naming an opaque class.  Catch it here and
    # say what to do.  (A list / numpy array / tuple of vectors has no .all_solutions/.get_tracker.)
    if hasattr(start_points, 'all_solutions') and hasattr(start_points, 'get_tracker'):
        raise TypeError(
            "HomotopySolver: start_points must be the actual start *points* (an iterable of "
            "solution vectors), but a {} solver was passed.  Call its .solve() and then pass "
            "its .all_solutions():\n"
            "    start_solver.solve()\n"
            "    nag_algorithm.HomotopySolver(homotopy, start_solver.all_solutions(), target)"
            .format(type(start_points).__name__))
    points = _coerce_start_points(start_points)
    # Validate the shapes HERE, at the Python boundary.  The native tracker checks them too, but
    # it throws from inside the tracking loop (possibly on a worker thread), where the exception
    # cannot propagate: the process aborts with no traceback, taking every computation in flight
    # with it (issues #369, #383).  A mismatch is a caller error, so it must be a Python exception.
    # The homotopy settles whether this solve is projective, because that was decided where it was
    # built (b2#382).  Bring the target into the SAME coordinates rather than making the caller do
    # it: homogenize a clone -- never their system, whose content is its identity and keys its
    # records -- and take the patch FROM the homotopy, since a fresh draw would be a different
    # patch and so different coordinates.
    if homotopy.is_homogeneous() and not target.is_homogeneous() and target.is_polynomial():
        from bertini import system as _bsys
        try:
            projective_target = _bsys.clone(target)
            projective_target.homogenize()
            projective_target.copy_patches(homotopy)
            target = projective_target
        except RuntimeError:
            pass        # cannot be written projectively; the consistency check below will speak

    n_vars = homotopy.num_variables()
    # A projective homotopy has more coordinates than the user writes -- one homogenizing
    # coordinate per affine group -- and its points live on a patch.  A caller holding solutions
    # from an affine solve has the shorter vectors, and lifting them is mechanical, so do it here
    # rather than make them find homogenize_point (b2#382).  A point that is already the full
    # length is left alone.
    import numpy as _np

    def _lift(p):
        if len(p) == n_vars:
            return p
        try:
            return homotopy.homogenize_point(_np.asarray(p))
        except Exception:
            return p        # not a liftable length either; the check below says so, in Python
    points = [_lift(p) for p in points]
    for k, p in enumerate(points):
        if len(p) != n_vars:
            raise ValueError(
                "HomotopySolver: start point {} has {} coordinate(s) but the homotopy has {} "
                "variable(s) (the path variable is not a coordinate); every start point needs "
                "exactly one coordinate per variable, in the homotopy's variable order"
                .format(k, len(p), n_vars))
    n_fun = homotopy.num_functions()
    if n_fun != n_vars:
        raise ValueError(
            "HomotopySolver: the homotopy has {} function(s) in {} variable(s), and tracking needs "
            "a SQUARE system.  Square it first -- randomize it down to the variable count "
            "(System.randomize) or cut it with a slice -- instead of tracking the overdetermined "
            "system.".format(n_fun, n_vars))
    # Adaptive precision derives its error bounds from the system's degree, and a homotopy that
    # is not polynomial does not have one.  The C++ side refuses too, but it would refuse while
    # the adaptive tracker is being constructed, and the text a user would see is about degree
    # bounds rather than about the choice they actually made.  Tracking analytic homotopies at
    # fixed precision works; only the adaptive criteria are unavailable.  See issue #439.
    if mptype == 'adaptive' and amp_config is None and not homotopy.is_polynomial():
        raise ValueError(
            "HomotopySolver: this homotopy is not polynomial, and adaptive precision needs a "
            "degree bound it therefore cannot have.  Tracking still works -- pass "
            "mptype='double' or mptype='multiple' to use a fixed-precision tracker.  If you want "
            "adaptive precision anyway, build a bertini.tracking.AMPConfig, set its "
            "jacobian_eval_error_bound and function_eval_error_bound to values you can defend for "
            "this system, and pass it as amp_config=.")

    user_start = _pybnalag.UserStartSystem(target, points)
    solver = solver_cls(target, user_start, homotopy)
    holder = _HomotopySolverHolder(solver, homotopy, target, user_start)

    if amp_config is not None and mptype == 'adaptive':
        holder.get_tracker().precision_setup(amp_config)

    return holder


def user_homotopy(homotopy, start_points, target, *, mptype='adaptive', precision=None, endgame='powerseries'):
    """Thin forwarder to :func:`HomotopySolver`, kept for back-compatibility."""
    return HomotopySolver(homotopy, start_points, target, mptype=mptype, precision=precision,
                          endgame=endgame)


def _projectivize_operands(systems):
    """Clones of the systems a homotopy is about to be built from, homogenized and patched.

    This is the only moment at which a homotopy can be made projective: once the blend exists its
    operands are fixed and cannot be homogenized, which is why the option lives on the builders
    rather than on the solver.  :func:`ZeroDimSolver` does the same thing at the same point, being
    a builder too.

    Clones, never the caller's systems.  A System's identity is its content, and records key a
    solve on that identity, so homogenizing one in place would silently change what the caller's
    own system *is* and stop an earlier solve of it being recalled.

    The patch is a random draw, so exactly one clone draws it and the rest copy: independent draws
    would put the operands in different coordinates.

    Parameters
    ----------
    systems : sequence of System
        The operand systems, sharing a variable structure.

    Returns
    -------
    list of System
        Projectivized clones, or the originals when this family cannot be written projectively.
    """
    # ADR-0020 and ADR-0026 fix a blend's operands at construction; b2#463 is the half-converted
    # state that used to result from trying to homogenize around one afterwards
    from bertini import system as _bsys

    live = [s for s in systems if s is not None]
    if not live:
        return list(systems)

    homogeneous = [s.is_homogeneous() for s in live]
    if any(homogeneous):
        # All of them: nothing to do.  SOME of them: leave it alone rather than quietly convert
        # the rest, and let the builder refuse.  A caller who hands over one projective system and
        # one affine one has made a mistake far more often than they have made a request, and
        # repairing it here would hide the mistake in the cases where it is one.
        return list(systems)
    if not all(s.is_polynomial() for s in live):
        return list(systems)    # homogenizing is about degrees, which these have not got

    # Not every block can homogenize -- a products-of-linears block over more than one affine group
    # refuses -- and a refusal partway through would leave some operands projective and some not.
    # Convert them all or hand back what we were given.
    clones = [_bsys.clone(s) for s in live]
    try:
        for s in clones:
            s.homogenize()
    except RuntimeError:
        return list(systems)    # this family cannot be written projectively; build it as authored

    clones[0].auto_patch()
    for s in clones[1:]:
        s.copy_patches(clones[0])

    out, it = [], iter(clones)
    for s in systems:
        out.append(next(it) if s is not None else None)
    return out


def straight_line_homotopy(target, start, *, path_variable='t', gamma=None, projectivize=True):
    """Form the straight-line homotopy H = (1-t)*target + gamma*t*start.

    The linear deformation from a start system to a target: at t=1 it is ``gamma*start``, whose
    solutions are your start points, and at t=0 it is ``target``.  Pair it with
    :func:`user_homotopy` or :class:`HomotopySolver` and the start system's solutions::

        start_solver = nag_algorithm.ZeroDimSolver(start, mptype='adaptive')
        start_solver.solve()
        H = nag_algorithm.straight_line_homotopy(target, start)
        solver = nag_algorithm.HomotopySolver(H, start_solver.all_solutions(), target)
        solver.solve()

    ``target`` and ``start`` must be built over the SAME variable objects, since the deformation
    combines their function trees.

    This also works when ``start`` carries a *structured evaluation block* -- a products-of-linears
    start built with :meth:`~bertini.System.add_products_of_linears`, say.  Such a block cannot be
    fused by System node arithmetic, which would silently drop it, so the two systems are combined
    with a blend block that evaluates whole Systems; that is the same construction the zero-dim
    solver uses for its own generated start systems (ADR-0020).

    Parameters
    ----------
    target : System
        The system whose solutions you want, reached at t=0.
    start : System
        The start system, whose known solutions are the start points.
    path_variable : str
        Name of the path variable t added to the homotopy (default ``'t'``).
    gamma : node, int or None
        The gamma coefficient.  ``None`` (default) draws a random one -- the gamma trick, which
        is what makes the path generic, so leave it alone unless you have a reason.  Pass an exact
        node (from :func:`bertini.coefficient`) off the real axis for a reproducible path.

        Pass ``1`` for a **coefficient-parameter homotopy**, where the deformation is a path in
        parameter space and no gamma belongs in it: then ``start`` is a generic member of the
        family, solved once and reused for any number of targets that share it.  Its coefficients
        should be genuinely generic (random complex), which is what keeps the straight path off
        the (measure-zero) singular locus.

    projectivize : bool, default True
        Homogenize and patch ``target`` and ``start`` before combining them, so the homotopy is
        tracked over projective coordinates and infinity is an ordinary place a path can reach
        instead of somewhere it runs off to.  This is what :func:`ZeroDimSolver` does with the
        systems it builds for itself.  It works on clones, so the systems you pass are untouched,
        and :func:`HomotopySolver` brings your affine target and start points into the homotopy's
        coordinates for you.  Pass False to build the homotopy affinely exactly as written -- the
        right choice when the affine coordinates are the point, as in all-real tracking.

        The two must be in the SAME coordinates to begin with: mixing a projective system with an
        affine one is refused rather than quietly repaired.

    Returns
    -------
    System
        The homotopy; pair it with :func:`user_homotopy` and your start points to solve.
    """
    from bertini._pybertini import system as _system
    if isinstance(gamma, int) and not isinstance(gamma, bool):
        from bertini.symbolics import Integer
        gamma = Integer(gamma)
    # a mismatched variable structure -- one end projective and the other affine, say -- is
    # refused by the native builder, which every caller goes through
    if projectivize:
        target, start = _projectivize_operands([target, start])
    return _system.make_homotopy(target, start, path_variable, gamma)


def moving_homotopy(fixed, start_moving, end_moving, *, path_variable='t', gamma=None, projectivize=True):
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

    Any of ``fixed``, ``start_moving`` and ``end_moving`` may be a :class:`~bertini.Slice`; it is
    taken as the system of its linear forms (a moving slice is the common case, issue #381).

    ``projectivize`` (default True) homogenizes and patches the three systems before they are
    combined, so the homotopy is tracked projectively and infinity is a place a path can reach --
    the same thing :func:`ZeroDimSolver` does with the systems it builds for itself.  It works on
    clones, so the systems you pass are untouched, and :func:`HomotopySolver` brings your affine
    target and start points into the homotopy's coordinates for you.  ``False`` builds the homotopy
    affinely exactly as written, which is what you want when the affine coordinates are the point,
    as in all-real tracking.

    The three must be in the SAME coordinates: mixing a projective system with an affine one is
    refused rather than quietly repaired, since it is a mistake far more often than a request.

    Returns the homotopy System; pair it with :func:`user_homotopy` and your start points to solve.
    """
    from bertini._pybertini import system as _system
    fixed, start_moving, end_moving = (s.as_system() if isinstance(s, _pybnalag.Slice) else s
                                       for s in (fixed, start_moving, end_moving))
    # as in straight_line_homotopy, the native builder refuses a mismatched variable structure
    if projectivize:
        fixed, start_moving, end_moving = _projectivize_operands([fixed, start_moving, end_moving])
    return _system.make_moving_homotopy(fixed, start_moving, end_moving, path_variable, gamma)


def parameter_sweep(make_system, generic_parameters, target_parameters,
                    *, collect=None, comm=None, mptype='adaptive', endgame='powerseries'):
    """Solve a whole family of systems that differ only in their coefficients.

    This is the *parameter homotopy* workhorse: pay for the hard ab-initio solve **once**, at a
    generic parameter value, then reach every parameter point you actually care about by cheap
    tracking that reuses those start solutions.  It is also the unit of parallelism -- the points
    are independent, so pass an MPI communicator and the sweep spreads across the ranks for you.

    Parameters
    ----------
    make_system : callable
        ``make_system(parameters) -> bertini.System``.  Must return systems of the SAME shape
        every call -- same variables and same monomials -- with only the coefficient *values*
        depending on ``parameters``.  (Until first-class coefficient parameters land, this factory
        is how you say "the same system at a different parameter value".)
    generic_parameters
        The parameter value for the one-time generic solve.  For robustness this should be
        *generic* -- random complex values -- so the straight-line coefficient path to each target
        avoids the measure-zero singular locus.
    target_parameters : iterable
        The parameter values you want solved.  Results line up with this order.
    collect : callable, optional
        ``collect(solver) -> value``.  If given, each solved point is reduced to ``collect(solver)``
        and the return is the list of those values (instead of the solvers).  **Required** when
        ``comm`` is given -- solver objects cannot cross MPI ranks, so what is gathered must be a
        picklable value (e.g. a solution count, or ``solver.to_dataframe()``).
    comm : mpi4py communicator, optional
        If given, ``target_parameters`` is split across the ranks (each rank solves the generic
        once, then tracks its slice -- with the path tracking inside each solve threaded across the
        rank's cores via ``OMP_NUM_THREADS``).  The collected results are gathered and **every rank
        returns the full list in the original order**.  This is the two-level model: MPI across
        parameter points, threads across paths -- one extra argument, no manual scatter/gather.
    mptype, endgame
        Passed through to the solves (see :func:`ZeroDim`).

    Returns
    -------
    list, one entry per ``target_parameters``: the solved :class:`ZeroDim` solvers, or -- if
    ``collect`` is given -- the ``collect(solver)`` values.  A solver exposes ``.all_solutions()``,
    ``.solution_metadata()`` and ``.to_dataframe()`` (let the solver do the real/finite
    classification -- don't redo it by hand).

    Examples
    --------
    Serial / threaded (threads are automatic -- a solve uses all cores by default)::

        solvers = parameter_sweep(make_system, generic, targets)
        counts  = [positive_real_count(s) for s in solvers]

    Across MPI ranks, threads within -- the user-facing MPI code is just ``comm=`` and ``collect=``::

        from mpi4py import MPI
        counts = parameter_sweep(make_system, generic, targets,
                                 collect=positive_real_count, comm=MPI.COMM_WORLD)
        # every rank now holds the full `counts` list, in `targets` order
    """
    targets = list(target_parameters)
    if comm is not None and collect is None:
        raise ValueError("parameter_sweep(comm=...) needs collect=<solver -> picklable value>: "
                         "solver objects cannot be gathered across MPI ranks.")

    if comm is None:
        my_indices = range(len(targets))
    else:
        my_indices = range(comm.Get_rank(), len(targets), comm.Get_size())

    generic = make_system(generic_parameters)
    gen_solver = ZeroDimSolver(generic, mptype=mptype, endgame=endgame)
    gen_solver.solve()
    start_points = gen_solver.all_solutions()

    local = []   # (index, solver-or-collected) for this rank's slice, in index order
    for i in my_indices:
        target = make_system(targets[i])
        # gamma=1: this is a coefficient-parameter homotopy, a path in parameter space, so the
        # gamma trick has no place in it
        H = straight_line_homotopy(target, generic, gamma=1)
        solver = HomotopySolver(H, start_points, target, mptype=mptype, endgame=endgame)
        solver.solve()
        local.append((i, collect(solver) if collect is not None else solver))

    if comm is None:
        return [value for _, value in local]

    # Gather every rank's (index, value) pairs and reassemble in the original target order, so
    # every rank returns the same full list.
    merged = {}
    for chunk in comm.allgather(local):
        merged.update(dict(chunk))
    return [merged[i] for i in range(len(targets))]


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
            # event.tracker() is the tracker that ACTUALLY runs this path: the solver's member
            # tracker in a serial solve, or the thread-local clone in a threaded solve.  Attaching
            # here (rather than to event.solver().get_tracker()) is what makes per-path collection
            # work identically with and without threads -- under threading the member tracker runs
            # nothing, so collecting from it would silently yield empty series.
            tracker = event.tracker()
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


# NID (numerical irreducible decomposition) is framework scaffolding whose Solve() is not yet
# implemented, and its witness-set datatypes ride with it -- hide them from the public surface
# (import * / tab-completion / docs) until it works.  They stay importable by explicit name.
_HIDDEN = tuple(n for n in dir(_pybnalag)
                if n.startswith(('NID', 'NumericalIrreducibleDecomposition', 'WitnessSet')))
__all__ = [n for n in dir(_pybnalag) if n not in _HIDDEN]
__all__.append('ZeroDimSolver')
__all__.append('HomotopySolver')
__all__.append('user_homotopy')
__all__.append('straight_line_homotopy')
__all__.append('parameter_sweep')
__all__.append('moving_homotopy')
__all__.append('SolutionPathCollector')


# DoublePrecisionTotalDegree = bertini._pybertini.nag_algorithms.ZeroDimCauchyDoublePrecisionTotalDegree
