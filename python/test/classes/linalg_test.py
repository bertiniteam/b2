"""The linear-algebra layer: vectors/matrices of variables and exact coefficients.

bertini.linalg lets you write conditions over vectors and matrices of variables with
ordinary numpy operators (A @ x - lam*x), with the hard rule that coefficients stay
exact -- Python floats are refused so they cannot poison the arbitrary-precision tree.
"""

from fractions import Fraction

import numpy as np
import pytest

import bertini as pb
from bertini import linalg as bla
from bertini import multiprec as mp


def test_variable_vector_names_and_shape():
    x = bla.variable_vector('x', 3)
    assert x.shape == (3,)
    assert [v.name for v in x] == ['x0', 'x1', 'x2']
    y = bla.variable_vector('y', 2, start=1)
    assert [v.name for v in y] == ['y1', 'y2']


def test_variable_matrix_names_and_shape():
    M = bla.variable_matrix('m', 2, 3)
    assert M.shape == (2, 3)
    assert M[0, 0].name == 'm_0_0'
    assert M[1, 2].name == 'm_1_2'


def test_coefficient_accepts_exact_values():
    # ints (python and numpy), Fractions, exact strings, and bertini multiprecision values
    for v in [2, np.int64(2), Fraction(3, 4), '2.5', '3/4',
              mp.Complex('1.5', '2.5'), mp.Float('1.5')]:
        node = bla.coefficient(v)
        # the result is a usable function-tree node: it combines with a variable
        _ = node * pb.Variable('z')
    # an existing node passes through unchanged
    z = pb.Variable('z')
    assert bla.coefficient(z) is z


@pytest.mark.parametrize("bad", [2.5, -0.1, np.float64(2.5), 1 + 2j, np.complex128(1j), True])
def test_coefficient_refuses_floats_and_bools(bad):
    # the whole point: a 64-bit float (or a bool masquerading as 1) must not silently
    # enter the arbitrary-precision tree.
    with pytest.raises(TypeError):
        bla.coefficient(bad)


def test_as_coefficients_matrix_of_exact_strings():
    A = bla.as_coefficients([['5/2', '1'], ['0', '3']])
    assert A.shape == (2, 2)
    x = bla.variable_vector('x', 2)
    eqs = A @ x                                  # builds two linear expressions
    assert eqs.shape == (2,)


def test_as_coefficients_refuses_a_float_anywhere():
    with pytest.raises(TypeError):
        bla.as_coefficients([[2, 1], [0, 3.0]])   # the lone 3.0 is rejected


def test_integer_matrix_times_variable_vector_needs_no_coercion():
    # numpy int * Variable already yields Integer coefficients, so an integer matrix
    # flows straight through @ with no as_coefficients call.
    A = np.array([[2, 1], [0, 3]])
    x = bla.variable_vector('x', 2)
    eqs = A @ x
    assert eqs.shape == (2,)


def test_add_functions_adds_each_component():
    x = bla.variable_vector('x', 3)
    A = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    sys = pb.System()
    n = bla.add_functions(sys, A @ x, basename='f')
    assert n == 3
    assert sys.num_functions() == 3


def test_add_functions_accepts_a_single_expression():
    x = bla.variable_vector('x', 2)
    sys = pb.System()
    assert bla.add_functions(sys, np.array([3, 4]) @ x - 1) == 1
    assert sys.num_functions() == 1
