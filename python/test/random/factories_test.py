"""Issue #294: friendly random factories (random_real / random_complex / random_vector) and the
``vg @ coeffs`` projection sugar.

A ``random_real`` is a genuine ``real_mp`` (not a ``complex_mp`` with zero imaginary part), the draws
are seed-reproducible and continuous (so a real projection direction is generic, unlike the quantized
orthonormal ``random_matrix``), and ``vg @ coeffs`` builds the single linear-combination node.
"""

import numpy as np

import bertini as pb
from bertini.multiprec import real_mp, complex_mp


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
