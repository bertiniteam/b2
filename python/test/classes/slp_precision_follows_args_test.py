"""A System must evaluate at whatever precision its arguments carry (bertini2 #377).

The compiled SLP is a precision-independent tape of operations; only the Memory holding
values carries digits.  So the Memory's precision is an artifact of the current evaluation,
never an invariant, and evaluation re-tags rather than refuses.

This used to throw, and there was no way out.  A Memory takes its precision from the
ambient ``default_precision()`` when the program is lazily compiled, while the owning
System keeps whatever it was told, so the two diverge the moment anything moves the ambient
default -- which an AMP tracker or endgame does as a matter of course.  Both setters
short-circuit when handed the value they already hold, so neither could repair it:

    RuntimeError: variable_values and SLP must be of same precision.
                  respective precisions: 30 16

Mirrors the C++ suite in core/test/classes/precision_propagation_test.cpp.
"""

import numpy as np
import pytest

import bertini
from bertini import System, VariableGroup
from bertini.symbolics import Rational, Variable
from bertini.multiprec import complex_mp, real_mp


def two_var_system():
    """f(x, y) = x^2 - y; Jacobian [2x, -1] is exact at any precision."""
    x, y = Variable("x"), Variable("y")
    s = System()
    vg = VariableGroup()
    vg.append(x)
    vg.append(y)
    s.add_variable_group(vg)
    s.add_functions([x**2 - y])
    return s


def point_at(prec, re):
    saved = bertini.default_precision()
    bertini.default_precision(prec)
    try:
        p = np.array([complex_mp(str(re), '0'), complex_mp('1', '0')])
    finally:
        bertini.default_precision(saved)
    return p


def test_eval_after_the_ambient_default_moved_underneath_the_system():
    """The wedge, reproduced: build at 30, let something move the default to 16, evaluate
    at the system's OWN precision.  Before the fix this raised, and no setter could undo
    it -- the System was unusable in both directions."""
    bertini.default_precision(30)
    s = two_var_system()

    bertini.default_precision(16)          # what an AMP tracker/endgame leaves behind

    v = np.atleast_1d(np.asarray(s.eval(point_at(30, 3))))
    assert abs(complex(float(repr(v[0].real)), float(repr(v[0].imag))) - 8) < 1e-25


def test_jacobian_after_the_ambient_default_moved_underneath_the_system():
    bertini.default_precision(30)
    s = two_var_system()
    bertini.default_precision(16)

    J = np.atleast_2d(np.asarray(s.eval_jacobian(point_at(30, 3))))
    assert J.shape == (1, 2)
    assert abs(float(repr(J[0, 0].real)) - 6.0) < 1e-25     # d/dx = 2x
    assert abs(float(repr(J[0, 1].real)) + 1.0) < 1e-25     # d/dy = -1


@pytest.mark.parametrize("ladder", [[16, 50, 30, 100, 20], [100, 16, 100]])
def test_one_system_evaluated_up_and_down_a_precision_ladder(ladder):
    """Precision is per-evaluation, so a single System serves any sequence of them."""
    bertini.default_precision(16)
    s = two_var_system()
    for prec in ladder:
        v = np.atleast_1d(np.asarray(s.eval(point_at(prec, 3))))
        assert abs(float(repr(v[0].real)) - 8.0) < 1e-14


def test_retag_rebuilds_constants_rather_than_zero_padding_them():
    """The mean one.  A constant that is INEXACT in decimal -- 1/3 -- must come back at
    full accuracy after the memory is re-tagged upward, which is only true if re-tagging
    refills constants from their exact recipes.  Padding a 16-digit 1/3 out to 100 digits
    leaves the tail zero, and the residual would be ~1e-17 instead of ~1e-100.
    """
    x = Variable("x")
    s = System()
    vg = VariableGroup()
    vg.append(x)
    s.add_variable_group(vg)
    # NOTE Rational(a, b) is (real, imaginary), not numerator/denominator -- the
    # single-argument spelling is the one that means one third
    s.add_functions([x - Rational('1/3')])

    bertini.default_precision(16)                 # force the lazy compile at 16 digits
    s.eval(np.array([complex_mp('1', '0')]))

    bertini.default_precision(100)
    third = np.array([complex_mp(real_mp(1) / real_mp(3), real_mp(0))])
    v = np.atleast_1d(np.asarray(s.eval(third)))
    residual = abs(float(repr(v[0].real)))
    assert residual < 1e-90, f"constant was padded, not rebuilt: residual {residual:g}"


def test_path_variable_precision_only_ever_increases_memory_precision():
    """The path variable arrives after the variables are already in memory, so re-tagging
    on it may only raise -- lowering would truncate the variables just written."""
    x, t = Variable("x"), Variable("t")
    s = System()
    vg = VariableGroup()
    vg.append(x)
    s.add_variable_group(vg)
    s.add_path_variable(t)
    s.add_functions([x * t])

    bertini.default_precision(30)
    pt = np.array([complex_mp('2', '0')])

    bertini.default_precision(60)
    time = complex_mp('3', '0')

    v = np.atleast_1d(np.asarray(s.eval(pt, time)))
    assert abs(float(repr(v[0].real)) - 6.0) < 1e-25
