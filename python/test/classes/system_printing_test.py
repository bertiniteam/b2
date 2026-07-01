"""Human-facing printing of block-composed systems.

print(system) / str(system) describes the system block by block.  A structured block keeps its
placeholder symbol for legibility (a linear form prints as ``c.[x, y, 1]``, a randomization as
``R . g``) and shows the actual coefficients just below it -- short (4 significant figures) in the
default/terse form, full working precision with system.describe(verbose=True).  Terse truncates a
block's listing after the first several rows so a large system does not flood the terminal.  The old
debugging leftovers -- the working-vector dump, the differentiated flag, the 'unnamed_function'
names -- are gone, and structured-block rows are no longer invisible.
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

    assert 'f_0 = x^2+y^2-1' in text
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
    assert 'R . g' in terse                       # placeholder label
    assert '(R: 2x3 randomization matrix)' in terse
    assert 'g_0 = x^2+y^2-1' in terse              # the underlying functions, indented and shown
    assert 'g_1 = x*y' in terse
    assert 'g_2 = x^2+y^2-x-y' in terse
    assert 'R =' in terse                          # the actual matrix is shown by default now
    assert terse.count('[') >= 2                   # its two rows

    verbose = r.describe(verbose=True)
    assert 'R =' in verbose                        # same layout, full precision
    for n in NOISE:
        assert n not in verbose


def test_linear_forms_block_placeholder_and_coefficients():
    x, y = pb.Variable('x'), pb.Variable('y')
    m = pb.System(); m.add_variable_group(_vg(x, y))
    m.add_function(x*x + y*y - 1)
    linalg.add_linear(m, np.array([[2, 1]]), np.array([x, y]), [-1])   # 2x + y - 1, a LinearFormsBlock

    terse = str(m)
    assert 'f_0 = x^2+y^2-1' in terse              # poly row
    assert 'f_1 = c.[x, y, 1]' in terse            # linear-forms placeholder kept (legible structure)
    assert 'c =' in terse                          # ... with the coefficients in a legend below
    assert '2' in terse and '-1' in terse          # the actual values (here exact integers)

    verbose = m.describe(verbose=True)
    assert 'c.[x, y, 1]' in verbose and 'c =' in verbose   # same layout, full precision
    assert '2' in verbose


def test_linear_form_coefficients_short_in_terse_full_in_verbose():
    import re
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    s = na.Slice.random_complex(_vg(x, y, z), 1).as_system()   # irrational random coefficients

    terse, verbose = s.describe(), s.describe(verbose=True)
    assert 'c.[x, y, z, 1]' in terse and 'c =' in terse

    # Count SIGNIFICANT FIGURES, not characters.  A character-width bound is magnitude-dependent
    # and flaky: a 4-significant-figure coefficient below 1e-3 (e.g. 0.0003162) is 9 characters
    # wide, so a `width <= 8` assertion fails on an unlucky small random draw.  Significant figures
    # are magnitude-independent -- terse is 4 sig figs for ANY coefficient -- so this never flakes.
    def most_sig_figs(text):
        return max((len(t.replace('.', '').lstrip('0')) for t in re.findall(r'\d+\.\d+', text)),
                   default=0)

    assert most_sig_figs(terse) <= 4               # terse: 4 significant figures, regardless of magnitude
    assert most_sig_figs(verbose) >= 12            # verbose: full working precision (~30 digits)


def test_terse_coefficient_width_varies_with_magnitude_but_sig_figs_do_not():
    # Deterministic companion to the test above (which uses random coefficients): with KNOWN
    # coefficients spanning magnitudes we pin the exact terse vs verbose rendering.  The small
    # coefficient (< 1e-3) renders as `0.0003162` -- 9 characters but still 4 significant figures.
    # This is the case that made the old `width <= 8` assertion flaky; it is correct behavior.
    import re
    pb.default_precision(30)
    x, y = pb.Variable('x'), pb.Variable('y')
    s = pb.System(); s.add_variable_group(_vg(x, y))
    s.add_function(x*x + y*y - 1)
    # exact decimal strings (add_linear rejects floats, which would cap block precision):
    # one coefficient below 1e-3, one of order 1, with many digits so verbose stays long.
    A = np.array([['0.00031622776601683794', '0.31622776601683794339']], dtype=object)
    linalg.add_linear(s, A, np.array([x, y]))

    terse, verbose = s.describe(), s.describe(verbose=True)

    # terse: 4 significant figures for BOTH magnitudes ...
    assert '0.0003162' in terse                    # small coeff: 9 characters wide, 4 sig figs
    assert '0.3162' in terse                        # order-1 coeff: 6 characters wide, 4 sig figs
    assert len('0.0003162') == 9                    # documents the width that broke `<= 8`
    sig = lambda t: len(t.replace('.', '').lstrip('0'))
    assert all(sig(t) <= 4 for t in re.findall(r'\d+\.\d+', terse))

    # verbose: the full exact strings survive (full working precision, not 4 sig figs)
    assert '0.00031622776601683794' in verbose
    assert '0.31622776601683794339' in verbose


def test_terse_truncates_many_forms_but_verbose_shows_all():
    x, y = pb.Variable('x'), pb.Variable('y')
    coeffs = [[i + 1, i + 2, i + 3] for i in range(12)]        # 12 affine forms on (x, y), > the cap of 10
    s = linalg.slice_from_coefficients(coeffs, [x, y]).as_system()

    terse, verbose = s.describe(), s.describe(verbose=True)
    assert terse.count('c.[') == 10                # capped at kTerseRowCap
    assert 'more form' in terse                    # with a truncation note
    assert verbose.count('c.[') == 12              # verbose shows every form
    assert 'more form' not in verbose


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
    assert 'f_0 = x^2+y^2-1' in terse              # the fixed row is a plain polynomial, shown
    assert '(1-t)*A' in terse and 'blend of 2 systems' in terse
    assert 'path variable: t' in terse
    assert 'f_1..f_1' not in terse                 # a single moving row reads 'f_1', not 'f_1..f_1'

    verbose = H.describe(verbose=True)
    assert 'A_0 = y-x' in verbose and 'B_0 = y' in verbose          # operand functions, indented
