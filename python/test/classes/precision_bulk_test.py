"""Bulk precision of a whole vector/matrix: bertini.precision(A) and bertini.precision(A, n).

Each mp number carries its own working precision.  ``bertini.precision`` reads or re-casts the
precision of an entire container in one call (delegating to the C++ Precision(container) /
Precision(container, digits)).  The setter is functional -- it returns a NEW array at the requested
precision and leaves the original untouched -- because eigenpy marshals mp arrays by copy.
"""

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import complex_mp, real_mp

# the autouse conftest fixture resets default precision to this baseline before each test
BASELINE = 30


def _complex_vec():
    return np.array([complex_mp('1'), complex_mp('2'), complex_mp('3')])


def _complex_mat():
    return np.array([[complex_mp('1'), complex_mp('2')],
                     [complex_mp('3'), complex_mp('4')]])


def test_get_precision_vector_and_matrix():
    assert pb.precision(_complex_vec()) == BASELINE
    assert pb.precision(_complex_mat()) == BASELINE
    assert pb.precision(np.array([real_mp('1'), real_mp('2')])) == BASELINE
    assert pb.precision(np.array([[real_mp('1'), real_mp('2')],
                                  [real_mp('3'), real_mp('4')]])) == BASELINE


def test_set_precision_is_functional_and_leaves_original_alone():
    v = _complex_vec()
    w = pb.precision(v, 80)
    assert pb.precision(w) == 80          # the copy is re-cast
    assert pb.precision(v) == BASELINE    # the original is untouched (functional form)
    assert w.dtype == np.dtype(complex_mp)
    assert w[0].precision == 80           # every element actually carries the new precision


def test_set_precision_matrix_keeps_shape():
    M = _complex_mat()
    N = pb.precision(M, 100)
    assert N.shape == (2, 2)
    assert pb.precision(N) == 100
    assert all(N[i, j].precision == 100 for i in range(2) for j in range(2))


def test_set_precision_real():
    v = np.array([real_mp('1'), real_mp('2')])
    w = pb.precision(v, 60)
    assert pb.precision(w) == 60
    assert w.dtype == np.dtype(real_mp)
    assert pb.precision(v) == BASELINE


def test_round_trip_precision():
    v = _complex_vec()
    up = pb.precision(v, 90)
    back = pb.precision(up, BASELINE)
    assert pb.precision(back) == BASELINE


def test_top_level_and_submodule_are_the_same():
    # bertini.precision is the hoisted bertini.multiprec.precision
    assert pb.precision is pb.multiprec.precision
