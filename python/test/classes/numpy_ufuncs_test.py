# This file is part of Bertini 2.
#
# python/test/classes/numpy_ufuncs_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/classes/numpy_ufuncs_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/classes/numpy_ufuncs_test.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   silviana amethyst
#   summer 2026

"""
Tests for the numpy ufunc coverage of the multiprecision dtypes.

Every registered loop calls the same boost::multiprecision free function that
the ``bertini.multiprec`` scalar function of the same name binds, so a ufunc
applied to an array must agree element-wise with the scalar function --
exactly, not approximately, since both run identical code at identical
precision.

Also pins:
- the identity-seeded reductions (np.sum/np.prod/np.mean) that older docs
  declared unsupported (they work on current numpy; regression-guard them),
- the sort/argsort/searchsorted/argmax family (dtype compare/arg slots),
- numpy semantics corners: rint half-to-even, remainder sign-of-divisor,
  nan propagation in minimum/maximum vs fmin/fmax,
- precision preservation through every loop shape,
- unwritten-slot safety (np.empty) per the ADR-0006 guard doctrine,
- the mp.imag copy-paste regression (returned the real part), and the array
  overloads of multiprec real/imag/arg that replace the silently-wrong
  ndarray .real/.imag attributes and np.angle.
"""

import numpy as np
import pytest

import bertini.multiprec as mp
from bertini.multiprec import complex_mp, real_mp


@pytest.fixture(params=[real_mp, complex_mp], ids=["real_mp", "complex_mp"])
def dtype(request):
    return request.param


def _sample(dtype):
    """A small array of values safe for every domain-restricted function."""
    if dtype is real_mp:
        return np.array([real_mp('0.25'), real_mp('0.5'), real_mp('0.75')])
    return np.array([complex_mp('0.25', '0.125'),
                     complex_mp('0.5', '-0.25'),
                     complex_mp('-0.75', '0.375')])


# ufuncs defined for both dtypes, with the matching multiprec scalar function
UFUNCS_BOTH = [
    (np.exp, mp.exp), (np.log, mp.log), (np.sqrt, mp.sqrt),
    (np.sin, mp.sin), (np.cos, mp.cos), (np.tan, mp.tan),
    (np.arcsin, mp.asin), (np.arccos, mp.acos), (np.arctan, mp.atan),
    (np.sinh, mp.sinh), (np.cosh, mp.cosh), (np.tanh, mp.tanh),
    (np.arcsinh, mp.asinh), (np.arccosh, None), (np.arctanh, mp.atanh),
]


class TestElementwiseAgreesWithScalarFunctions:
    """np.f(arr)[i] must equal mp.f(arr[i]) exactly (identical code path)."""

    @pytest.mark.parametrize(
        "ufunc,scalar", UFUNCS_BOTH,
        ids=[u.__name__ for u, _ in UFUNCS_BOTH])
    def test_transcendental(self, dtype, ufunc, scalar):
        if ufunc is np.arccosh:
            # acosh needs |x| >= 1 on the real line
            v = (np.array([real_mp('1.5'), real_mp(2)]) if dtype is real_mp
                 else np.array([complex_mp('1.5', '0.5')]))
            scalar = mp.acosh
        else:
            v = _sample(dtype)
        out = ufunc(v)
        assert out.dtype == np.dtype(dtype)
        for got, x in zip(out, v):
            assert got == scalar(x)

    def test_absolute(self, dtype):
        v = _sample(dtype)
        out = np.abs(v)
        # output is ALWAYS real, also for complex input (the magnitude)
        assert out.dtype == np.dtype(real_mp)
        for got, x in zip(out, v):
            assert got == mp.abs(x)

    def test_conjugate(self, dtype):
        v = _sample(dtype)
        out = np.conj(v)
        assert out.dtype == np.dtype(dtype)
        for got, x in zip(out, v):
            assert got == (mp.conj(x) if dtype is complex_mp else x)

    def test_power(self, dtype):
        v = _sample(dtype)
        out = np.power(v, v)
        for got, x in zip(out, v):
            assert got == x ** x

    def test_reciprocal(self, dtype):
        v = _sample(dtype)
        out = np.reciprocal(v)
        for got, x in zip(out, v):
            assert got == dtype(1) / x

    def test_square_negative_positive(self, dtype):
        v = _sample(dtype)
        assert all(np.square(v)[i] == v[i] * v[i] for i in range(len(v)))
        assert all(np.negative(v)[i] == -v[i] for i in range(len(v)))
        assert all(np.positive(v)[i] == v[i] for i in range(len(v)))

    def test_real_only_transcendentals(self):
        v = _sample(real_mp)
        pairs = [(np.log10, mp.log), (np.exp2, None), (np.log2, None),
                 (np.expm1, None), (np.log1p, None), (np.cbrt, None)]
        # spot-check values against the double versions loosely; the exact
        # contract (same boost call as a scalar) has no bound scalar twin for
        # these, so compare against float64 at double precision.
        for ufunc, _ in pairs:
            got = ufunc(v)
            want = ufunc(np.array([float(x) for x in v]))
            for g, w in zip(got, want):
                assert abs(float(g) - w) < 1e-14, ufunc.__name__

    def test_sign_real(self):
        v = np.array([real_mp(-3), real_mp(0), real_mp(2)])
        assert [str(s) for s in np.sign(v)] == ['-1', '0', '1']

    def test_sign_complex_is_unit_modulus(self):
        # numpy-2 semantics: sign(z) = z/|z|, 0 at 0
        z = complex_mp(3, 4)
        s = np.sign(np.array([z, complex_mp(0)]))
        assert s[0] == complex_mp('0.6', '0.8')
        assert s[1] == complex_mp(0)

    def test_arctan2_hypot_copysign(self):
        # atan2/hypot use dedicated algorithms, so agreement with the composed
        # formulas is to the last ulp, not bit-exact
        tol = real_mp('1e-25')
        a = np.array([real_mp(1), real_mp(-2)])
        b = np.array([real_mp(3), real_mp(4)])
        assert mp.abs(np.arctan2(a, b)[0] - mp.atan(a[0] / b[0])) < tol
        assert mp.abs(np.hypot(a, b)[1] - mp.sqrt(a[1] * a[1] + b[1] * b[1])) < tol
        got = np.copysign(a, np.array([real_mp(-1), real_mp(1)]))
        assert [str(x) for x in got] == ['-1', '2']


class TestNumpySemanticsCorners:
    """Corners where numpy's semantics differ from naive C/boost calls."""

    def test_rint_rounds_half_to_even(self):
        # regression: boost's rint rounds half AWAY from zero; numpy (and the
        # loop, via mpfr_rint in MPFR_RNDN) rounds half to even
        v = np.array([real_mp('0.5'), real_mp('1.5'), real_mp('2.5'),
                      real_mp('-0.5'), real_mp('-2.5')])
        assert [str(x) for x in np.rint(v)] == ['0', '2', '2', '-0', '-2']

    def test_rint_and_round_on_complex(self):
        # numpy defines rint (and hence np.round) for complex, component-wise
        w = np.array([complex_mp('1.5', '2.5'), complex_mp('-0.5', '3.4')])
        want = np.rint(np.array([1.5 + 2.5j, -0.5 + 3.4j]))
        for got, ref in zip(np.rint(w), want):
            assert complex(got) == ref
        for got, ref in zip(np.round(w), want):
            assert complex(got) == ref

    def test_floor_ceil_trunc(self):
        v = np.array([real_mp('1.7'), real_mp('-1.7')])
        assert [str(x) for x in np.floor(v)] == ['1', '-2']
        assert [str(x) for x in np.ceil(v)] == ['2', '-1']
        assert [str(x) for x in np.trunc(v)] == ['1', '-1']

    def test_remainder_takes_sign_of_divisor(self):
        a = np.array([real_mp(7), real_mp(-7), real_mp(7), real_mp(-7)])
        b = np.array([real_mp(3), real_mp(3), real_mp(-3), real_mp(-3)])
        assert [str(x) for x in np.mod(a, b)] == ['1', '2', '-2', '-1']
        # fmod keeps C semantics (sign of dividend)
        assert [str(x) for x in np.fmod(a, b)] == ['1', '-1', '1', '-1']
        assert [str(x) for x in np.floor_divide(a, b)] == ['2', '-3', '-3', '2']

    def test_minimum_maximum_propagate_nan_fmin_fmax_ignore_it(self):
        nan = real_mp('nan')
        one = np.array([real_mp(1)])
        nans = np.array([nan])
        assert mp.abs(np.fmax(nans, one)[0] - real_mp(1)) == 0
        assert mp.abs(np.fmin(nans, one)[0] - real_mp(1)) == 0
        assert np.isnan(np.maximum(nans, one))[0]
        assert np.isnan(np.minimum(nans, one))[0]

    def test_minimum_maximum_values(self):
        a = np.array([real_mp(1), real_mp(5)])
        b = np.array([real_mp(3), real_mp(2)])
        assert [str(x) for x in np.minimum(a, b)] == ['1', '2']
        assert [str(x) for x in np.maximum(a, b)] == ['3', '5']

    def test_predicates(self, dtype):
        good = _sample(dtype)
        assert not np.isnan(good).any()
        assert np.isfinite(good).all()
        assert not np.isinf(good).any()
        assert np.isnan(good).dtype == np.dtype(bool)
        if dtype is real_mp:
            bad = np.array([real_mp('nan'), real_mp('inf'), real_mp(1)])
        else:
            bad = np.array([complex_mp('nan', '0'), complex_mp('0', 'inf'),
                            complex_mp(1, 1)])
        assert list(np.isnan(bad)) == [True, False, False]
        assert list(np.isinf(bad)) == [False, True, False]
        assert list(np.isfinite(bad)) == [False, False, True]

    def test_signbit(self):
        v = np.array([real_mp(-2), real_mp(0), real_mp(3)])
        assert list(np.signbit(v)) == [True, False, False]


class TestSortingAndArgExtrema:
    """The dtype compare/argmax/argmin slots (real only -- complex is unordered)."""

    def test_sort_and_argsort(self):
        v = np.array([real_mp(3), real_mp(1), real_mp(2)])
        assert [str(x) for x in np.sort(v)] == ['1', '2', '3']
        assert list(np.argsort(v)) == [1, 2, 0]

    def test_argmax_argmin_max_min(self):
        v = np.array([real_mp(3), real_mp(1), real_mp(7), real_mp(2)])
        assert np.argmax(v) == 2
        assert np.argmin(v) == 1
        assert str(np.max(v)) == '7'
        assert str(np.min(v)) == '1'

    def test_identityless_reduce_does_not_corrupt_input(self):
        # numpy >= 2.5 initializes the accumulator of an identityless reduce
        # (np.min/np.max) by BITWISE copy of element 0, so the accumulator and
        # v[0] share one mpfr allocation.  The loops must never free or write
        # through an output slot's existing allocation (slot_write): before
        # that rule, np.max silently rewrote v[0] and freed its storage, and
        # the next reduce crashed the interpreter (use-after-free -> corrupted
        # allocator).  Values AND the input array must survive, repeatedly.
        v = np.array([real_mp(3), real_mp(1), real_mp(7), real_mp(2)])
        for _ in range(3):
            assert str(np.max(v)) == '7'
            assert [str(x) for x in v] == ['3', '1', '7', '2']
            assert str(np.min(v)) == '1'
            assert [str(x) for x in v] == ['3', '1', '7', '2']

    def test_argmax_nan_wins(self):
        # numpy float semantics: the first nan is the arg-extremum
        v = np.array([real_mp(1), real_mp('nan'), real_mp(3)])
        assert np.argmax(v) == 1
        assert np.argmin(v) == 1

    def test_searchsorted_and_median(self):
        v = np.array([real_mp(1), real_mp(2), real_mp(4)])
        assert np.searchsorted(v, real_mp(3)) == 2
        assert str(np.median(v)) == '2'

    def test_complex_stays_unordered(self):
        w = np.array([complex_mp(1, 2), complex_mp(0, 1)])
        with pytest.raises(TypeError):
            np.sort(w)


class TestReductions:
    """Identity-seeded reductions -- previously documented as crashing.

    They work on current numpy; these tests exist so any numpy/eigenpy
    combination that breaks them again fails loudly here instead of in
    user code.
    """

    def test_sum_prod_mean_real(self):
        v = np.array([real_mp(1), real_mp(2), real_mp(3)])
        assert np.sum(v) == real_mp(6)
        assert np.prod(v) == real_mp(6)
        assert np.mean(v) == real_mp(2)

    def test_sum_prod_mean_complex(self):
        w = np.array([complex_mp(1, 2), complex_mp(3, 4)])
        assert np.sum(w) == complex_mp(4, 6)
        assert np.prod(w) == complex_mp(-5, 10)
        assert np.mean(w) == complex_mp(2, 3)

    def test_bare_reduce_and_cumsum(self, dtype):
        v = np.array([dtype(1), dtype(2), dtype(3)])
        assert np.add.reduce(v) == dtype(6)
        assert list(np.cumsum(v)) == [dtype(1), dtype(3), dtype(6)]

    def test_reduce_with_explicit_initial_still_works(self, dtype):
        # the old explicit-initial workaround must keep working (portable to old numpy)
        v = np.array([dtype(1), dtype(2)])
        assert np.add.reduce(v, initial=dtype(10)) == dtype(13)


class TestCloseness:
    """The float64 boundary: tolerance ORDERINGS against a double are allowed
    (a comparison yields a bool -- no float flows into a multiprecision value;
    this mirrors the C++ solvers' double ToleranceT and the scalar
    GreatLessVisitor<T, double>).  Everything that would let a float VALUE into
    an mp computation stays closed: mixed equality, mixed arithmetic,
    np.isclose's internal float tolerances."""

    def test_tolerance_comparison_with_float(self, dtype):
        v = _sample(dtype)
        w = v + dtype('1e-20')
        assert np.all(np.abs(v - w) <= 1e-10)
        assert not np.all(np.abs(v - (w + dtype(1))) <= 1e-10)
        # both operand orders, and float64 arrays as well as scalars
        assert np.all(1e-10 >= np.abs(v - w))
        assert np.all(np.abs(v - w) < np.full(len(v), 1e-10))

    def test_tolerance_comparison_is_exact_not_sloppy(self):
        # the double is compared exactly (boost mixed compare), not by rounding
        # the mp value down to double first
        tiny = real_mp('1e-22')
        assert np.all(np.array([tiny]) < 1e-10)
        assert not np.any(np.array([tiny]) < 1e-30)

    def test_allclose_idiom_all_mp_still_works(self, dtype):
        v = _sample(dtype)
        w = v + dtype('1e-20')
        assert np.all(np.abs(v - w) <= real_mp('1e-10'))

    def test_float_equality_stays_blocked(self, dtype):
        # exact equality against a float literal is the 0.1-intent trap; it is
        # not bound at the scalar level either
        v = _sample(dtype)
        with pytest.raises(TypeError):
            v == 0.1

    def test_float_arithmetic_stays_blocked(self, dtype):
        v = _sample(dtype)
        with pytest.raises(TypeError):
            v + 0.1

    def test_isclose_itself_still_raises(self, dtype):
        # its internal float64 rtol/atol cannot promote; if this ever starts
        # passing, numpy grew user-dtype promotion -- revisit the numpy docs page
        v = _sample(dtype)
        with pytest.raises(TypeError):
            np.isclose(v, v)


class TestExplicitDownConversion:
    """astype is the explicit, conscious truncation to double precision."""

    def test_astype_float_and_complex(self):
        v = np.array([real_mp('1.5'), real_mp(2)])
        w = np.array([complex_mp(1, 2), complex_mp(3, 4)])
        assert list(v.astype(float)) == [1.5, 2.0]
        assert list(v.astype(complex)) == [1.5 + 0j, 2.0 + 0j]
        assert list(w.astype(complex)) == [1 + 2j, 3 + 4j]

    def test_astype_int_stays_forbidden(self, dtype):
        v = np.array([dtype(1)])
        with pytest.raises(TypeError):
            v.astype(np.int64)


class TestPrecisionPreservation:
    """Outputs carry the operands' precision, not the ambient default --
    including through the mixed real/complex division inside sign and
    reciprocal (the known boost precision-mis-tagging hazard)."""

    HIGH = 50

    def _high_precision_sample(self, dtype):
        mp.default_precision(self.HIGH)
        v = (np.array([real_mp('1.5')]) if dtype is real_mp
             else np.array([complex_mp('1.5', '2.5')]))
        mp.default_precision(30)
        assert v[0].precision == self.HIGH
        return v

    @pytest.mark.parametrize("ufunc", [np.exp, np.sqrt, np.abs, np.sign,
                                       np.reciprocal, np.conj, np.rint],
                             ids=lambda u: u.__name__)
    def test_unary_output_precision(self, dtype, ufunc):
        if dtype is complex_mp and ufunc is np.rint:
            pytest.skip("rint is real-only")
        v = self._high_precision_sample(dtype)
        assert ufunc(v)[0].precision == self.HIGH

    def test_binary_output_precision(self, dtype):
        v = self._high_precision_sample(dtype)
        assert np.power(v, v)[0].precision == self.HIGH


class TestUnwrittenSlotSafety:
    """Every new loop shape must survive never-written np.empty slots
    (which hold the all-zero BMP sentinel) -- the ADR-0006 doctrine."""

    def test_unary_loops_on_empty(self, dtype):
        e = np.empty(3, dtype=dtype)
        for ufunc in (np.exp, np.sin, np.conj, np.sign, np.abs,
                      np.isnan, np.isfinite):
            ufunc(e)  # must not crash

    def test_binary_loops_on_empty(self, dtype):
        e = np.empty(3, dtype=dtype)
        np.power(e, e)
        if dtype is real_mp:
            np.minimum(e, e)
            np.arctan2(e, e)
            np.mod(e, e)

    def test_sort_and_argmax_on_empty(self):
        e = np.empty(4, dtype=real_mp)
        np.sort(e)
        np.argmax(e)
        np.argmin(e)

    def test_reductions_on_zeros(self, dtype):
        z = np.zeros(3, dtype=dtype)
        assert np.sum(z) == dtype(0)


class TestScalarsAreOwnedCopies:
    """Regression tests for the getitem aliasing fix.

    getitem used to return boost::ref into the numpy buffer (as stock eigenpy
    does), so a scalar extracted from a temporary array -- most visibly the
    result of np.sum/np.mean -- dangled once the array was freed: reading it
    later gave zeros/garbage or SIGABRT inside mpfr (the ADR-0031 / #259
    hazard class).  getitem now returns an owned copy.
    """

    def test_reduce_scalar_survives_its_array(self, dtype):
        s = np.sum(np.array([dtype(1), dtype(2), dtype(3)]))
        # the source (temporary) array is gone; s must still be intact
        assert s == dtype(6)
        assert str(s) is not None  # printing used to MPFR-assert on the corpse
        assert s / dtype(3) == dtype(2)

    def test_mean_is_correct_not_silently_zero(self, dtype):
        # np.mean's internal divide ran on a dangling extraction and returned 0
        v = np.array([dtype(1), dtype(2), dtype(3)])
        assert np.mean(v) == dtype(2)

    def test_stored_indexed_elements_stay_distinct(self):
        # the ADR-0031 shape, now safe at the binding level (still copy in
        # Python code by convention)
        pts = np.array([complex_mp(1, 1), complex_mp(2, 2), complex_mp(3, 3)])
        kept = [pts[i] for i in range(3)]
        del pts
        assert [str(k) for k in kept] == ['(1,1)', '(2,2)', '(3,3)']

    def test_mutating_an_extracted_scalar_leaves_the_array_alone(self):
        v = np.array([real_mp(1), real_mp(2)])
        x = v[0]
        x += real_mp(10)
        assert v[0] == real_mp(1)


class TestComponentAccessors:
    """The multiprec real/imag/arg array overloads, and the scalar imag
    regression."""

    def test_scalar_imag_returns_imaginary_part(self):
        # regression: mp.imag was bound to boost::multiprecision::real by a
        # copy-paste error, so it returned the REAL part
        z = complex_mp(1, 2)
        assert mp.real(z) == real_mp(1)
        assert mp.imag(z) == real_mp(2)

    def test_array_real_imag(self):
        w = np.array([complex_mp(1, 2), complex_mp(3, 4)])
        r, i = mp.real(w), mp.imag(w)
        assert r.dtype == np.dtype(real_mp) and i.dtype == np.dtype(real_mp)
        assert [str(x) for x in r] == ['1', '3']
        assert [str(x) for x in i] == ['2', '4']

    def test_array_arg_replaces_np_angle(self):
        w = np.array([complex_mp(1, 1), complex_mp(-1, 0)])
        a = mp.arg(w)
        assert a.dtype == np.dtype(real_mp)
        assert a[0] == mp.arg(w[0])
        assert a[1] == mp.arg(w[1])

    def test_ndarray_real_imag_attributes_are_untrustworthy(self):
        # documenting-by-test: numpy cannot know a legacy user dtype is
        # complex-like, so ndarray .real returns the complex values themselves
        # and .imag returns zeros.  If numpy ever fixes this, the accessors
        # above become optional and the numpy docs page should be updated.
        w = np.array([complex_mp(1, 2)])
        assert w.real.dtype == np.dtype(complex_mp)   # not real_mp!
        assert w.imag[0] == complex_mp(0)             # wrong value, by numpy


class TestGuardedNumpyComponentFunctions:
    """np.real/np.imag/np.angle raise on plain mp-complex arrays instead of
    silently returning wrong values (bertini._numpy_guard) -- a crash is better
    than incorrect values.  Solution points override .real/.imag at the subclass
    level and pass through correct."""

    def test_np_real_imag_raise_on_plain_complex_mp_array(self):
        w = np.array([complex_mp(1, 2), complex_mp(3, 4)])
        with pytest.raises(TypeError, match="bertini.real"):
            np.real(w)
        with pytest.raises(TypeError, match="bertini.real"):
            np.imag(w)
        with pytest.raises(TypeError, match="bertini.real"):
            np.angle(w)

    def test_np_real_imag_raise_on_lists_of_complex_mp(self):
        # a list converts to a plain mp array inside numpy, same wrong path
        with pytest.raises(TypeError):
            np.imag([complex_mp(1, 2)])

    def test_guard_passes_everything_else_through(self):
        # ordinary numpy is untouched
        z = np.array([1 + 2j, 3 + 4j])
        assert list(np.real(z)) == [1.0, 3.0]
        assert list(np.imag(z)) == [2.0, 4.0]
        assert np.angle(np.array([1j]))[0] == pytest.approx(np.pi / 2)
        # real_mp arrays are not complex: base semantics are already correct
        v = np.array([real_mp(1), real_mp(2)])
        assert list(np.real(v)) == [real_mp(1), real_mp(2)]
        assert list(np.imag(v)) == [real_mp(0), real_mp(0)]
        # mp-complex SCALARS go through the (correct) scalar properties
        assert np.real(complex_mp(1, 2)) == real_mp(1)
        assert np.imag(complex_mp(1, 2)) == real_mp(2)

    def test_guard_is_idempotent(self):
        import bertini._numpy_guard as guard
        before = np.real
        guard.install()
        assert np.real is before

    def test_solution_real_imag_are_correct(self):
        from bertini.records import Solution
        s = Solution(np.array([complex_mp(1, 2), complex_mp(3, 4)]))
        assert [str(x) for x in s.real] == ['1', '3']
        assert [str(x) for x in s.imag] == ['2', '4']
        assert s.real.dtype == np.dtype(real_mp)
        # np.real/np.imag on a Solution route through the subclass property
        assert [str(x) for x in np.real(s)] == ['1', '3']
        assert [str(x) for x in np.imag(s)] == ['2', '4']

    def test_solution_real_imag_correct_for_double_solves_too(self):
        from bertini.records import Solution
        s = Solution(np.array([1 + 2j, 3 + 4j]))
        assert list(s.real) == [1.0, 3.0]
        assert list(s.imag) == [2.0, 4.0]

    def test_np_angle_raises_helpfully_for_all_mp_complex(self):
        # np.angle branches on the DTYPE (never the .real/.imag attributes), so
        # not even the Solution subclass can make it work -- the guard turns the
        # cryptic arctan2 failure into a pointer at mp.arg, for every spelling
        from bertini.records import Solution
        for val in (np.array([complex_mp(1, 1)]),
                    Solution(np.array([complex_mp(1, 1)])),
                    complex_mp(1, 1)):
            with pytest.raises(TypeError, match="arg"):
                np.angle(val)
