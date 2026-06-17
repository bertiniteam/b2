"""Human-facing printing of block-composed systems.

print(system) / str(system) describes the system block by block, with placeholder symbols for the
structured blocks' coefficients (terse, the default) or the actual numbers and underlying functions
(system.describe(verbose=True)).  The old debugging leftovers -- the working-vector dump, the
differentiated flag, the 'unnamed_function' names -- are gone, and structured-block rows are no
longer invisible.
"""

import numpy as np
import pytest

import bertini as pb
from bertini import linalg, nag_algorithm as na

NOISE = ('current variable values', 'not differentiated', 'is differentiated', 'unnamed_function')


def _vg(*vs):
    return pb.VariableGroup(list(vs))


def test_plain_polynomial_is_clean():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y))
    s.add_function(x*x + y*y - 1); s.add_function(x - y)
    text = str(s)

    assert 'f_0 = x*x+y*y-1' in text
    assert 'f_1 = x-y' in text
    for n in NOISE:
        assert n not in text          # debugging leftovers are gone
    assert text.count('= ') == 2      # exactly two functions printed (no extras, none missing)


def test_randomization_shows_placeholder_and_underlying_functions():
    x, y = pb.Variable('x'), pb.Variable('y')
    o = pb.System(); o.add_variable_group(_vg(x, y))
    o.add_function(x*x + y*y - 1); o.add_function(x*y); o.add_function(x*x + y*y - x - y)
    r = linalg.randomize(o)

    terse = str(r)
    assert 'R . g' in terse                       # placeholder, not the matrix
    assert '(R: 2x3 randomization matrix)' in terse
    assert 'g_0 = x*x+y*y-1' in terse              # the underlying functions, indented and shown
    assert 'g_1 = x*y' in terse
    assert 'g_2 = x*x+y*y-x-y' in terse
    assert 'R =' not in terse                      # the actual matrix is NOT in the terse form

    verbose = r.describe(verbose=True)
    assert 'R =' in verbose                        # ... but it is in verbose
    assert verbose.count('[') >= 2                 # the two matrix rows
    for n in NOISE:
        assert n not in verbose


def test_linear_forms_block_placeholder_vs_actual():
    x, y = pb.Variable('x'), pb.Variable('y')
    m = pb.System(); m.add_variable_group(_vg(x, y))
    m.add_function(x*x + y*y - 1)
    linalg.add_linear(m, np.array([[2, 1]]), np.array([x, y]), [-1])   # 2x + y - 1, a LinearFormsBlock

    terse = str(m)
    assert 'f_0 = x*x+y*y-1' in terse              # poly row
    assert 'f_1 = c.[x, y, 1]' in terse            # linear-forms row: placeholder, both rows visible

    verbose = m.describe(verbose=True)
    assert '(2)*x' in verbose and '(1)*y' in verbose and '(-1)' in verbose   # actual coefficients


def test_products_of_linears_block():
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y))
    linalg.add_products_of_linears(s, [[[1, 0, -1], [1, 0, 1]]])    # (x-1)(x+1)

    assert 'prod of 2 linear forms' in str(s)                       # terse placeholder
    verbose = s.describe(verbose=True)
    assert '*' in verbose and 'x' in verbose                        # actual product of factors


def test_moving_homotopy_blend_and_path_variable():
    x, y = pb.Variable('x'), pb.Variable('y')
    fx = pb.System(); fx.add_variable_group(_vg(x, y)); fx.add_function(x*x + y*y - 1)
    sm = pb.System(); sm.add_variable_group(_vg(x, y)); sm.add_function(y)
    em = pb.System(); em.add_variable_group(_vg(x, y)); em.add_function(y - x)
    H = na.moving_homotopy(fx, sm, em, gamma=linalg.coefficient(pb.multiprec.Complex('0.6', '0.8')))

    terse = str(H)
    assert 'f_0 = x*x+y*y-1' in terse              # the fixed row is a plain polynomial, shown
    assert '(1-t)*A' in terse and 'blend of 2 systems' in terse
    assert 'path variable: t' in terse
    assert 'f_1..f_1' not in terse                 # a single moving row reads 'f_1', not 'f_1..f_1'

    verbose = H.describe(verbose=True)
    assert 'A_0 = y-x' in verbose and 'B_0 = y' in verbose          # operand functions, indented
