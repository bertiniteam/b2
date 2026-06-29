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
    >>> from bertini import linalg                                       # doctest: +SKIP
    >>> import bertini                                                   # doctest: +SKIP
    >>> x, y = bertini.Variable('x'), bertini.Variable('y')             # doctest: +SKIP
    >>> linalg.slice_from_coefficients([[2, 3, 1]], [x, y]).coefficients().shape  # doctest: +SKIP
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
    >>> from bertini import linalg                                       # doctest: +SKIP
    >>> import bertini                                                   # doctest: +SKIP
    >>> x, y = bertini.Variable('x'), bertini.Variable('y')             # doctest: +SKIP
    >>> s = linalg.slice_from_coefficients([[2, 3, 1], [1, -1, 4]], [x, y])  # doctest: +SKIP
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
``bertini.linalg.slice_from_coefficients`` for exact numpy/list coefficients).
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
    'is_finite', 'is_real', 'is_singular', 'is_nonsolution', 'multiplicity', 'multiplicity_representative',
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
    """Give every bound solver class (ZeroDimSolver*, HomotopySolver*) a to_dataframe method (idempotent)."""
    for name in dir(_pybnalag):
        if not name.startswith(('ZeroDimSolver', 'HomotopySolver')):
            continue
        cls = getattr(_pybnalag, name)
        if isinstance(cls, type) and not getattr(cls, '_b2_has_to_dataframe', False):
            cls.to_dataframe = _zerodim_to_dataframe
            cls._b2_has_to_dataframe = True


_attach_to_dataframe()


# --- ZeroDimSolver: a friendly factory over the bound ZeroDimSolver<endgame x precision> classes ---

# Each bound solver class is named ZeroDimSolver<Endgame><Precision> (the start system is NO
# longer part of the type -- it is chosen at construction).  Select endgame + precision by string to
# pick the class, then pick the start system separately (default total_degree, else a factory).
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
    'totaldegree': 'total_degree', 'total_degree': 'total_degree', 'td': 'total_degree',
    'rootsofunity': 'roots_of_unity', 'roots_of_unity': 'roots_of_unity', 'rou': 'roots_of_unity',
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
    else -- two or more variable groups, or any homogeneous group -- is multihomogeneous (total
    degree / roots of unity would be the wrong, over-counting start system there).

    INTERIM: the single-affine-group case infers ``rootsofunity`` (the safe default).  The eventual
    inference is the linear-product ``totaldegree`` (general position), gated on the Cauchy-endgame
    fix; see core/.../nag_algorithms/common/policies.hpp.
    """
    if system.num_variable_groups() == 1 and system.num_hom_variable_groups() == 0:
        return 'rootsofunity'
    return 'mhom'


def ZeroDimSolver(system, *, endgame='cauchy', mptype='adaptive', startsystem='infer',
                  precision=None):
    """Construct a zero-dim solver by name, with friendly defaults.

    ``ZeroDimSolver(system)`` is the Cauchy endgame in adaptive precision with the start system
    **inferred from the system's variable-group structure** -- total degree for a single affine group,
    multihomogeneous otherwise -- so a multi-group (e.g. eigenvalue) system gets MHom automatically
    rather than an over-counting total-degree start.  Override any piece with a string::

        ZeroDimSolver(system, endgame='cauchy', mptype='amp', startsystem='mhom')

    Parameters
    ----------
    system : the polynomial :class:`~bertini.System` to solve.
    endgame : ``'cauchy'`` (default) or ``'powerseries'``.
    mptype : the precision -- ``'double'``, ``'multiple'``, or ``'adaptive'`` (``'amp'``, the default).
    precision : an alias for ``mptype``; if given (not ``None``) it overrides ``mptype``.
    startsystem : ``'infer'`` (default -- choose from the variable-group structure, matching the
        C++ blackbox), or force it with ``'totaldegree'`` / ``'mhom'``.  To run from a homotopy you
        built yourself with given start points, use :class:`HomotopySolver` / :func:`blend_homotopy`
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
    if precision is not None:
        mptype = precision
    start_key = str(startsystem).strip().lower().replace('-', '_').replace('_', '')
    # user-homotopy can't be built from a system alone -- point at the right entry point.
    if start_key in ('user', 'userhomotopy'):
        raise ValueError(
            "ZeroDimSolver does not build the homotopy solver (its construction needs a homotopy "
            "and start points, not just a system); build the homotopy with "
            "nag_algorithm.blend_homotopy / coefficient_parameter_homotopy and solve it with "
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
    # total_degree is the default constructor; any other start is selected with a factory.
    if start_enum == _pybnalag.StartSystemType.total_degree:
        return cls(system)
    return cls(system, _pybnalag.start_system_factory(start_enum))


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


def HomotopySolver(homotopy, start_points, target, *, precision='adaptive', endgame='cauchy'):
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
        when the variable coordinates line up (e.g. an affine homotopy).
    target : System
        The system the solutions satisfy at t=0 -- used for dehomogenize / residual and for the
        solver's consistency check.  It must NOT have a path variable.
    precision : {'adaptive', 'double', 'multiple'}
        'adaptive' (default) is the robust path.
    endgame : {'cauchy', 'powerseries'}

    Returns a solver: call ``.solve()`` then ``.all_solutions()`` as for any solver.
    """
    prec = {'double': 'double', 'multiple': 'multiple', 'fixed_multiple': 'multiple',
            'adaptive': 'adaptive'}.get(precision, precision)
    eg = {'cauchy': 'cauchy', 'powerseries': 'powerseries',
          'power_series': 'powerseries'}.get(endgame, endgame)
    try:
        cls_name = _HOMOTOPY_SOLVER_CLASSES[(prec, eg)]
    except KeyError:
        raise ValueError(
            "HomotopySolver: unknown (precision, endgame) = ({!r}, {!r}); "
            "precision in {{'adaptive','double','multiple'}}, endgame in {{'cauchy','powerseries'}}"
            .format(precision, endgame))
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
    user_start = _pybnalag.UserStartSystem(target, list(start_points))
    solver = solver_cls(target, user_start, homotopy)
    return _HomotopySolverHolder(solver, homotopy, target, user_start)


def user_homotopy(homotopy, start_points, target, *, precision='adaptive', endgame='cauchy'):
    """Thin forwarder to :func:`HomotopySolver`, kept for back-compatibility."""
    return HomotopySolver(homotopy, start_points, target, precision=precision, endgame=endgame)


def coefficient_parameter_homotopy(target, generic, path_variable='t'):
    """Build a parameter homotopy interpolating two systems of the same shape.

    Returns H = (1 - t) * target + t * generic with ``t`` added as its path variable, so at t=1
    it is ``generic`` (whose solutions are your start points) and at t=0 it is ``target``.  Pair
    it with :func:`user_homotopy`: solve ``generic`` once, then reuse its solutions to move to
    ``target`` (and to any number of further targets that share ``generic``)::

        gen_solver = nag_algorithm.ZeroDimSolver(generic, mptype='adaptive')
        gen_solver.solve()
        H = nag_algorithm.coefficient_parameter_homotopy(target, generic)
        solver = nag_algorithm.HomotopySolver(H, gen_solver.all_solutions(), target)
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


def parameter_sweep(make_system, generic_parameters, target_parameters,
                    *, collect=None, comm=None, mptype='adaptive', endgame='cauchy'):
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
        H = coefficient_parameter_homotopy(target, generic)
        solver = HomotopySolver(H, start_points, target, precision=mptype, endgame=endgame)
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


__all__ = dir(_pybnalag)
__all__.append('ZeroDimSolver')
__all__.append('HomotopySolver')
__all__.append('user_homotopy')
__all__.append('coefficient_parameter_homotopy')
__all__.append('parameter_sweep')
__all__.append('moving_homotopy')
__all__.append('blend_homotopy')
__all__.append('SolutionPathCollector')


# DoublePrecisionTotalDegree = bertini._pybertini.nag_algorithms.ZeroDimCauchyDoublePrecisionTotalDegree