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

# config structs gain update()/repr/to_dict/...; algorithm classes gain configure().
from ..config import enhance_all, enhance_owners
enhance_all(_pybnalag)
enhance_owners(_pybnalag)


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


def ZeroDim(system, *, endgame='cauchy', mptype='multiple', startsystem='totaldegree',
            precision=None):
    """Construct a zero-dim solver by name, with friendly defaults.

    ``ZeroDim(system)`` is the Cauchy endgame, multiple precision, total-degree start system --
    i.e. ``ZeroDimCauchyFixedMultiplePrecisionTotalDegree(system)`` -- without typing that out.
    Select the other 17 combinations with strings::

        ZeroDim(system, endgame='cauchy', mptype='amp', startsystem='mhom')

    Parameters
    ----------
    system : the polynomial :class:`~bertini.System` to solve.
    endgame : ``'cauchy'`` (default) or ``'powerseries'``.
    mptype : the precision -- ``'double'``, ``'multiple'`` (default), or ``'adaptive'`` (``'amp'``).
    precision : an alias for ``mptype``; if given (not ``None``) it overrides ``mptype``.
    startsystem : ``'totaldegree'`` (default) or ``'mhom'``.  To run from a homotopy you built
        yourself with given start points, use :func:`user_homotopy` / :func:`blend_homotopy`
        instead (their construction needs the homotopy and start points, not just a system).

    Returns a solver; call ``.solve()`` then ``.solutions()`` as for any zero-dim solver.
    """
    if precision is not None:
        mptype = precision
    # user-homotopy can't be built from a system alone -- point at the right entry point.
    if str(startsystem).strip().lower().replace('-', '_').replace('_', '') in ('user', 'userhomotopy'):
        raise ValueError(
            "ZeroDim does not build the user-homotopy solver (its construction needs a homotopy "
            "and start points, not just a system); build the homotopy with "
            "nag_algorithm.blend_homotopy / coefficient_parameter_homotopy and solve it with "
            "nag_algorithm.user_homotopy(homotopy, start_points, target).")
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

    Returns a solver: call ``.solve()`` then ``.solutions()`` as for any zero-dim solver.
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
        solver = nag_algorithm.user_homotopy(H, gen_solver.solutions(), target)
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

    Returns the homotopy System; pair it with :func:`user_homotopy` and your start points to solve.
    """
    from bertini._pybertini import system as _system
    return _system.make_homotopy(target, start, path_variable, gamma)


__all__ = dir(_pybnalag)
__all__.append('ZeroDim')
__all__.append('user_homotopy')
__all__.append('coefficient_parameter_homotopy')
__all__.append('blend_homotopy')


# DoublePrecisionTotalDegree = bertini._pybertini.nag_algorithms.ZeroDimCauchyDoublePrecisionTotalDegree