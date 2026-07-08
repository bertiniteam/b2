"""Vectorized mp helpers (#298, #301): real/imag/abs/conj/round/sum/norm/is_real stay mp-native over
a scalar / list / numpy array, so values do not collapse to float64."""

import numpy as np
import pytest

import bertini as pb
from bertini.multiprec import real_mp, complex_mp


@pytest.fixture
def pt():
    return np.array([complex_mp("1.5"), complex_mp(2, 3)], dtype=object)


def test_real_imag_stay_real_mp(pt):
    re = pb.real(pt)
    im = pb.imag(pt)
    assert all(isinstance(v, real_mp) for v in re)
    assert all(isinstance(v, real_mp) for v in im)
    assert [str(v) for v in re] == ['1.5', '2']
    assert [str(v) for v in im] == ['0', '3']


def test_abs_is_real_mp_magnitude(pt):
    a = pb.abs(pt)
    assert all(isinstance(v, real_mp) for v in a)
    assert abs(float(a[1]) - (13 ** 0.5)) < 1e-12   # |2 + 3i| = sqrt(13)


def test_conj_negates_imaginary(pt):
    c = pb.conj(pt)
    assert complex(c[1]) == complex(2, -3)


def test_round_stays_mp(pt):
    r = pb.round(real_mp("2.34567"), 2)
    assert isinstance(r, real_mp)
    assert str(r) == '2.35'
    # complex rounds both parts
    rc = pb.round(complex_mp("1.23456", "7.89123"), 2)
    assert isinstance(rc, complex_mp)


def test_is_real_predicate(pt):
    assert pb.is_real(pt) is False                                    # has a 2+3i entry
    assert pb.is_real([complex_mp("1.0"), complex_mp("2.0")]) is True
    assert pb.is_real([complex_mp(0, 1e-20)], tol=1e-10) is True      # tiny imaginary within tol


def test_sum_and_norm_stay_mp():
    s = pb.sum([real_mp(1), real_mp(2), real_mp(3)])
    assert isinstance(s, real_mp)
    assert float(s) == 6.0
    n = pb.norm([real_mp(3), real_mp(4)])
    assert isinstance(n, real_mp)
    assert abs(float(n) - 5.0) < 1e-40                                # full-precision 5, not float noise


def test_helpers_work_on_scalars_too():
    assert isinstance(pb.real(complex_mp(1, 2)), real_mp)
    assert isinstance(pb.abs(complex_mp(3, 4)), real_mp)
    assert abs(float(pb.abs(complex_mp(3, 4))) - 5.0) < 1e-12
