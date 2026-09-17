"""str is the readable Python spelling, repr the exact one, to_classic() Bertini 1's (ADR-0059).

Interface checks: repr round-trips through eval in the bertini namespace, str never shows a
caret or a (re,im) pair, and the classic spelling is what to_classic_input writes.  The
exact spellings themselves are pinned in C++ (print_dialect_test.cpp).
"""

import numpy as np
import pytest

import bertini
from bertini import System, Variable
from bertini.symbolics import Rational, Integer, Complex
from bertini.multiprec import real_mp, complex_mp


@pytest.fixture
def xy():
    return Variable('x'), Variable('y')


def _namespace(**names):
    ns = dict(vars(bertini))
    ns.update(vars(bertini.symbolics))      # Complex, Rational, Integer, ... by their node names
    ns.update(names)
    return ns


def _values_agree_exactly(a, b, x, y):
    s, t = System(), System()
    for sys_, f in ((s, a), (t, b)):
        sys_.add_variable_group([x, y])
        sys_.add_function(f)
    pt = np.array([complex_mp('0.37', '-1.25'), complex_mp('2.5', '0.125')])
    va = np.asarray(s.eval(pt)).ravel()[0]
    vb = np.asarray(t.eval(pt)).ravel()[0]
    return va == vb


def test_str_is_python_and_to_classic_is_bertini1(xy):
    x, y = xy
    e = x**2 + Rational('1/3', '1/2') * y - 3
    assert str(e) == 'x**2+(1/3+1/2*I)*y-3'
    assert e.to_classic() == 'x^2+(1/3+1/2*I)*y-3'
    assert '^' not in str(e) and '**' not in e.to_classic()


@pytest.mark.parametrize("precision", [30, 50], indirect=True)
def test_repr_rebuilds_an_equal_node_at_full_precision(xy, precision):
    x, y = xy
    third = real_mp(1) / real_mp(3)
    c = Complex(complex_mp(third, -third))              # a genuinely complex mp constant
    e = 3 * x**2 + Rational('1/3') * y - c * x + bertini.coefficient(third) * y**3 + Integer(7)

    text = repr(e)
    assert '**' in text and '^' not in text
    assert f"', {precision})" in text                   # the stored precision travels in the spelling

    rebuilt = eval(text, _namespace(x=x, y=y))
    assert repr(rebuilt) == text                        # the spelling is a fixed point
    assert _values_agree_exactly(rebuilt, e, x, y)      # and the values are bit-identical


def test_repr_of_a_bare_multiprecision_constant_is_exact_at_a_different_default():
    bertini.default_precision(50)
    c = bertini.coefficient(real_mp('0.1'))            # stored at 50 digits
    text = repr(c)
    assert text.startswith("real_mp('") and text.endswith("', 50)")
    bertini.default_precision(16)                       # a coarser session default ...
    assert complex_mp(eval(text, _namespace())) == c.value()     # ... still rebuilds the same value


def test_system_str_is_python_and_classic_input_is_bertini1(xy):
    x, y = xy
    s = System()
    s.add_variable_group([x, y])
    s.add_function(x**2 + y**2 - 1)
    s.add_function(Rational('1/3', '1/2') * x - y)
    assert 'x**2+y**2-1' in str(s)
    assert '(1/3+1/2*I)*x-y' in str(s)
    classic = s.to_classic_input()
    assert 'x^2+y^2-1' in classic and '(1/3+1/2*I)*x-y' in classic
    # the classic spelling reads straight back (to the same VALUES: Bertini 1 has no exact
    # rationals, so 1/3 comes back as a division at the working precision)
    back = bertini.parse.system(classic)
    pt = np.array([complex(0.37, -1.25), complex(2.5, 0.125)])
    assert np.allclose(np.asarray(back.eval(pt)), np.asarray(s.eval(pt)), rtol=1e-14, atol=0)
