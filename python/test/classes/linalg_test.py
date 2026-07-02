"""The linear-algebra layer: vectors/matrices of variables and exact coefficients.

The System.add_* methods, together with bertini.variables / bertini.coefficient(s), let you
write conditions over vectors and matrices of variables with ordinary numpy operators
(A @ x - lam*x), with the hard rule that coefficients stay
exact -- Python floats are refused so they cannot poison the arbitrary-precision tree.
"""

from fractions import Fraction

import numpy as np
import pytest

import bertini as pb
from bertini import multiprec as mp


def test_variable_vector_names_and_shape():
    x = np.array(pb.variables('x', 3), dtype=object)
    assert x.shape == (3,)
    assert [v.name for v in x] == ['x0', 'x1', 'x2']
    y = np.array(pb.variables('y', range(1, 3)), dtype=object)
    assert [v.name for v in y] == ['y1', 'y2']


def test_variable_matrix_names_and_shape():
    M = np.array([[pb.Variable(f'm_{i}_{j}') for j in range(3)] for i in range(2)], dtype=object)
    assert M.shape == (2, 3)
    assert M[0, 0].name == 'm_0_0'
    assert M[1, 2].name == 'm_1_2'


def test_coefficient_accepts_exact_values():
    # ints (python and numpy), Fractions, exact strings, and bertini multiprecision values
    for v in [2, np.int64(2), Fraction(3, 4), '2.5', '3/4',
              mp.Complex('1.5', '2.5'), mp.Float('1.5')]:
        node = pb.coefficient(v)
        # the result is a usable function-tree node: it combines with a variable
        _ = node * pb.Variable('z')
    # an existing node passes through unchanged
    z = pb.Variable('z')
    assert pb.coefficient(z) is z


@pytest.mark.parametrize("bad", [2.5, -0.1, np.float64(2.5), 1 + 2j, np.complex128(1j), True])
def test_coefficient_refuses_floats_and_bools(bad):
    # the whole point: a 64-bit float (or a bool masquerading as 1) must not silently
    # enter the arbitrary-precision tree.
    with pytest.raises(TypeError):
        pb.coefficient(bad)


def test_as_coefficients_matrix_of_exact_strings():
    A = pb.coefficients([['5/2', '1'], ['0', '3']])
    assert A.shape == (2, 2)
    x = np.array(pb.variables('x', 2), dtype=object)
    eqs = A @ x                                  # builds two linear expressions
    assert eqs.shape == (2,)


def test_as_coefficients_refuses_a_float_anywhere():
    with pytest.raises(TypeError):
        pb.coefficients([[2, 1], [0, 3.0]])   # the lone 3.0 is rejected


def test_integer_matrix_times_variable_vector_needs_no_coercion():
    # numpy int * Variable already yields Integer coefficients, so an integer matrix
    # flows straight through @ with no as_coefficients call.
    A = np.array([[2, 1], [0, 3]])
    x = np.array(pb.variables('x', 2), dtype=object)
    eqs = A @ x
    assert eqs.shape == (2,)


def test_add_functions_adds_each_component():
    x = np.array(pb.variables('x', 3), dtype=object)
    A = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    sys = pb.System()
    n = sys.add_functions(A @ x)
    assert n == 3
    assert sys.num_functions() == 3


def test_add_functions_accepts_a_single_expression():
    x = np.array(pb.variables('x', 2), dtype=object)
    sys = pb.System()
    assert sys.add_functions(np.array([3, 4]) @ x - 1) == 1
    assert sys.num_functions() == 1


# --- the C++ LinearFormsBlock, driven from Python via add_linear_forms ---

def test_add_linear_forms_block_evaluates():
    # f0 = 2x + 3y + 1,  f1 = x - y + 4  (augmented rows; trailing column is the constant)
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_linear_forms([[2, 3, 1], [1, -1, 4]])
    assert sys.num_functions() == 2

    v = sys.eval(np.array([mp.Complex('1'), mp.Complex('1')], dtype=mp.Complex))
    assert abs(complex(v[0]) - 6) < 1e-10
    assert abs(complex(v[1]) - 4) < 1e-10


def test_add_linear_forms_accepts_exact_nonintegers():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    # 5/2 x + 1 y + 0   evaluated at (2, 1) = 5 + 1 = 6
    sys.add_linear_forms([['5/2', '1', '0']])
    v = sys.eval(np.array([mp.Complex('2'), mp.Complex('1')], dtype=mp.Complex))
    assert abs(complex(v[0]) - 6) < 1e-10


def test_add_linear_forms_refuses_floats():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    with pytest.raises(TypeError):
        sys.add_linear_forms([[2.5, 1, 0]])


# --- add_linear: the auto-target (A @ x + b = 0 as a LinearFormsBlock) ---

def test_add_linear_matches_scalar_expansion():
    # Build the SAME linear system two ways on systems with identical variable structure,
    # and check they evaluate identically: scalar function-tree expansion vs LinearFormsBlock.
    A = [[2, 3], [1, -1], [0, 4]]
    b = [1, 4, -2]

    # scalar version: add_functions(A @ x + b)
    xs = np.array(pb.variables('x', 2), dtype=object)
    sys_scalar = pb.System()
    sys_scalar.add_variable_group(pb.VariableGroup(list(xs)))
    sys_scalar.add_functions(np.array(A) @ xs + np.array(b))

    # block version: add_linear(A, x, b)
    xb = np.array(pb.variables('x', 2), dtype=object)
    sys_block = pb.System()
    sys_block.add_variable_group(pb.VariableGroup(list(xb)))
    sys_block.add_linear(A, xb, b)

    pt = np.array([mp.Complex('2'), mp.Complex('-1')], dtype=mp.Complex)
    v_scalar = sys_scalar.eval(pt)
    v_block = sys_block.eval(pt)

    assert len(v_scalar) == len(v_block) == 3
    for i in range(3):
        assert abs(complex(v_scalar[i]) - complex(v_block[i])) < 1e-10


def test_add_linear_refuses_floats():
    xs = np.array(pb.variables('x', 2), dtype=object)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(xs)))
    with pytest.raises(TypeError):
        sys.add_linear([[2.5, 1]], xs)


# --- products of linears: the first-class C++ ProductsOfLinearsBlock from Python ---

def test_add_products_of_linears_evaluates_and_degrees():
    # f0 = (x + 1)(x - 1) = x^2 - 1   (two linear factors -> degree 2)
    # f1 = 2x + 3y + 1                (a single factor      -> degree 1)
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_products_of_linears([
        [[1, 0, 1], [1, 0, -1]],   # rows are factors; the trailing column is the constant
        [[2, 3, 1]],
    ])
    assert sys.num_functions() == 2
    # a product's degree is its number of factors -- the point of the products-of-linears block.
    assert list(sys.degrees()) == [2, 1]

    v = sys.eval(np.array([mp.Complex('2'), mp.Complex('1')], dtype=mp.Complex))
    assert abs(complex(v[0]) - 3) < 1e-10   # (2 + 1)(2 - 1) = 3
    assert abs(complex(v[1]) - 8) < 1e-10   # 2*2 + 3*1 + 1 = 8


def test_add_products_of_linears_accepts_exact_nonintegers():
    # one function, one factor (1/2)x + 3/4; at x = 1 -> 1/2 + 3/4 = 5/4
    x = pb.Variable('x')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x]))
    sys.add_products_of_linears([[['1/2', '3/4']]])
    v = sys.eval(np.array([mp.Complex('1')], dtype=mp.Complex))
    assert abs(complex(v[0]) - 1.25) < 1e-10


def test_add_products_of_linears_refuses_floats():
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    with pytest.raises(TypeError):
        sys.add_products_of_linears([[[2.5, 0, 1], [1, 0, -1]]])


def test_add_products_of_linears_survives_clone():
    # System pickling isn't exposed to Python, but bertini.system.clone is the deep-copy
    # round-trip (the C++ boost-archive round-trip is covered by system_blocks_test).  A
    # products-of-linears block must survive it: coefficients copied, and the copy independent.
    x, y = pb.Variable('x'), pb.Variable('y')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup([x, y]))
    sys.add_products_of_linears([
        [[1, 0, 1], [1, 0, -1]],   # (x + 1)(x - 1) = x^2 - 1
        [[2, 3, 1]],               # 2x + 3y + 1
    ])

    clone = pb.system.clone(sys)
    assert list(clone.degrees()) == [2, 1]

    pt = np.array([mp.Complex('2'), mp.Complex('1')], dtype=mp.Complex)
    v0, v1 = sys.eval(pt), clone.eval(pt)
    assert len(v0) == len(v1) == 2
    for a, b in zip(v0, v1):
        assert abs(complex(a) - complex(b)) < 1e-12

    # independence: changing the original after cloning does not touch the clone.
    sys.add_function(x + y)
    assert sys.num_functions() == 3
    assert clone.num_functions() == 2


def test_add_products_of_linears_projective_homogeneous_forms():
    # Authoring works for a projective (homogeneous) variable group too: use homogeneous linear
    # factors (constant column zero).  (x0 - x1)(x0 + x1) -> 0 at (1,1), 8 at (3,1).
    x0, x1 = pb.Variable('x0'), pb.Variable('x1')
    sys = pb.System()
    sys.add_hom_variable_group(pb.VariableGroup([x0, x1]))
    sys.add_products_of_linears([[[1, -1, 0], [1, 1, 0]]])
    assert list(sys.degrees()) == [2]
    assert abs(complex(sys.eval(np.array([mp.Complex('1'), mp.Complex('1')], dtype=mp.Complex))[0])) < 1e-12
    assert abs(complex(sys.eval(np.array([mp.Complex('3'), mp.Complex('1')], dtype=mp.Complex))[0]) - 8) < 1e-10


# --- degrees of the linear-algebra evaluation paths ---------------------------------------
# After the polynomial-path fold, System.degrees() asks the blocks.  These pin the degree of
# each linear-algebra construction, with the eigenvalue distinction (bilinear vs linear) front
# and centre.  list(...) also confirms the returned std::vector<int> converts cleanly.

def test_degrees_add_linear_is_one():
    # A @ x = 0 as a LinearFormsBlock -- each row is degree 1.
    x = np.array(pb.variables('x', 2), dtype=object)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_linear(np.array([[2, 1], [0, 3]]), x)
    assert list(sys.degrees()) == [1, 1]


def test_degrees_add_linear_forms_is_one():
    # the augmented-matrix entry point -- also a LinearFormsBlock, degree 1.
    x = np.array(pb.variables('x', 2), dtype=object)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_linear_forms([[2, 1, -1], [0, 3, 4]])   # rows are (num_vars + 1) wide
    assert list(sys.degrees()) == [1, 1]


def test_degrees_scalar_linear_via_add_functions_is_one():
    # A @ x - 1 expanded to scalar function-tree rows is still degree 1 (polynomial block).
    x = np.array(pb.variables('x', 2), dtype=object)
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_functions(np.array([[2, 1], [0, 3]]) @ x - 1)
    assert list(sys.degrees()) == [1, 1]


def test_degrees_eigenvalue_bilinear_is_two():
    # (A - lambda I) x has a lambda*x term -> each row is degree 2.  This is the distinction
    # that makes the eigenvalue problem genuinely nonlinear (vs a constant-coefficient solve).
    x = np.array(pb.variables('x', 2), dtype=object)
    lam = pb.Variable('lam')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_variable_group(pb.VariableGroup([lam]))
    sys.add_functions(np.array([[2, 1], [1, 3]]) @ x - lam * x)
    assert list(sys.degrees()) == [2, 2]


def test_degrees_eigenvalue_with_normalization_is_mixed():
    # the full eigenvalue formulation: two bilinear rows (degree 2) plus a constant-coefficient
    # normalization carried as a LinearFormsBlock (degree 1) -> {2, 2, 1}, in block order.
    x = np.array(pb.variables('x', 2), dtype=object)
    lam = pb.Variable('lam')
    sys = pb.System()
    sys.add_variable_group(pb.VariableGroup(list(x)))
    sys.add_variable_group(pb.VariableGroup([lam]))
    sys.add_functions(np.array([[2, 1], [1, 3]]) @ x - lam * x)   # degree-2 rows
    sys.add_linear(np.array([[1, 1]]), x)                         # degree-1 normalization
    assert list(sys.degrees()) == [2, 2, 1]
