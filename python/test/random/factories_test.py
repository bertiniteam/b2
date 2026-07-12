"""Issue #294: friendly random factories (random_real / random_complex / random_vector) and the
``vg @ coeffs`` projection sugar.

A ``random_real`` is a genuine ``real_mp`` (not a ``complex_mp`` with zero imaginary part), the draws
are seed-reproducible and continuous (so a real projection direction is generic, unlike the quantized
orthonormal ``random_matrix``), and ``vg @ coeffs`` builds the single linear-combination node.
"""

import numpy as np

import bertini as pb
from bertini import symbolics as sym
from bertini.multiprec import real_mp, complex_mp
from bertini.symbolics import AbstractNode, Complex


def test_random_real_is_real_mp():
    r = pb.random_real()
    assert isinstance(r, real_mp)


def test_random_complex_is_complex_mp():
    c = pb.random_complex()
    assert isinstance(c, complex_mp)


def test_random_vector_element_types():
    vr = pb.random_vector(4, real=True)
    vc = pb.random_vector(4, real=False)
    assert len(vr) == 4 and len(vc) == 4
    assert all(isinstance(x, real_mp) for x in vr)
    assert all(isinstance(x, complex_mp) for x in vc)


def test_random_vector_seed_reproducible_and_seed_sensitive():
    pb.random.set_random_seed(1234)
    a = [repr(x) for x in pb.random_vector(5, real=True)]

    pb.random.set_random_seed(1234)
    b = [repr(x) for x in pb.random_vector(5, real=True)]
    assert a == b                                  # same seed -> same vector

    pb.random.set_random_seed(4321)
    c = [repr(x) for x in pb.random_vector(5, real=True)]
    assert a != c                                  # a different seed -> a different vector


def test_variablegroup_matmul_is_the_linear_combination():
    x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')
    vg = pb.VariableGroup([x, y, z])
    coeffs = pb.random_vector(3, real=True)

    pi = vg @ coeffs                               # sum_i coeffs[i] * vg[i], a single node

    pt = [complex_mp(2), complex_mp(-3), complex_mp(5)]
    got = complex(pi.eval(x=pt[0], y=pt[1], z=pt[2]))
    want = complex(coeffs[0]) * 2 + complex(coeffs[1]) * (-3) + complex(coeffs[2]) * 5
    assert abs(got - want) < 1e-9, (got, want)


# --- symbolic (node) forms: random constants ready to drop into an expression ---

def test_random_real_symbolic_is_a_real_node():
    n = pb.random_real(symbolic=True)
    assert isinstance(n, AbstractNode)
    assert isinstance(n, Complex)                  # a real literal is a zero-imaginary Complex leaf
    assert n.value().imag == 0                      # ... and its imaginary part really is 0


def test_random_complex_symbolic_is_a_complex_node():
    n = pb.random_complex(symbolic=True)
    assert isinstance(n, AbstractNode)
    assert isinstance(n, Complex)


def test_symbolic_default_is_still_the_numeric_value():
    # symbolic defaults False -- the pre-existing value-returning behavior is unchanged
    assert isinstance(pb.random_real(), real_mp)
    assert isinstance(pb.random_complex(), complex_mp)
    assert isinstance(pb.random_real(symbolic=False), real_mp)


def test_random_vector_symbolic_is_an_array_of_nodes():
    v = pb.random_vector(4, symbolic=True)
    assert len(v) == 4
    assert all(isinstance(e, AbstractNode) for e in v)


def test_symbolics_namespace_shortcuts_return_nodes():
    # bertini.symbolics.random_real / random_complex are the always-symbolic spellings
    assert isinstance(sym.random_real(), Complex)
    assert isinstance(sym.random_complex(), Complex)
    assert sym.random_real().value().imag == 0


def test_symbolic_node_evaluates_to_the_value_it_wraps():
    # the node is a genuine constant carrying the drawn value: it evaluates back to it
    pb.random.set_random_seed(99)
    n = pb.random_real(symbolic=True)
    drawn = n.value()                               # the stored complex_mp
    got = n.eval()                                  # a constant node needs no point
    assert abs(complex(got) - complex(drawn)) < 1e-12


def test_symbolic_node_drops_into_an_expression():
    x = pb.Variable('x')
    a = pb.random_real(symbolic=True)
    expr = a * x                                    # a random real coefficient times a variable
    # eval at x = 4 must equal (coefficient) * 4
    got = complex(expr.eval(x=complex_mp(4)))
    want = complex(a.value()) * 4
    assert abs(got - want) < 1e-12, (got, want)


def test_symbolic_precision_follows_default_precision():
    # a random float node is only as precise as its draw: the current default precision
    pb.default_precision(200)
    assert pb.random_real(symbolic=True).value().precision == 200
    assert pb.random_complex(symbolic=True).value().precision == 200
    assert sym.random_real().value().precision == 200


def test_symbolic_draws_are_seed_reproducible():
    pb.random.set_random_seed(2024)
    a = repr(pb.random_real(symbolic=True))
    pb.random.set_random_seed(2024)
    b = repr(pb.random_real(symbolic=True))
    assert a == b                                   # same seed -> same node
    pb.random.set_random_seed(2025)
    c = repr(pb.random_real(symbolic=True))
    assert a != c                                   # a different seed -> a different node


def test_exact_rational_random_alternative_is_a_rational_node():
    # the precision-independent alternative advertised in the docstrings: an exact random rational
    from bertini.symbolics import Rational
    assert isinstance(Rational.rand_real(), Rational)
    assert isinstance(Rational.rand(), Rational)
