//This file is part of Bertini 2.
//
//test/flint/flint_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/flint/flint_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/flint/flint_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// What b2 relies on from FLINT (bertiniteam/b2#412), checked on every platform CI builds:
//   - ball arithmetic is rigorous where b2 will use it (products, quotients, division by a
//     ball around zero), and mpfr <-> arf conversion is exact;
//   - FLINT and b2 share one GMP and one MPFR;
//   - zero-filled memory is an exact zero (what numpy hands a dtype that needs initialization,
//     so the ADR-0006 crash cannot happen for a FLINT-backed type);
//   - a complex rectangle of Boost.Interval<double> encloses Arb's evaluation of the same
//     polynomial, the double interval type's design.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Bertini 2 FLINT"

#include <boost/test/unit_test.hpp>

// gmp.h and mpfr.h must come first: FLINT declares its mpz and mpfr conversions only when
// those headers are already in.
#include <gmp.h>
#include <mpfr.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arf.h>
#include <flint/arb.h>
#include <flint/acb.h>

#include <boost/numeric/interval.hpp>

#include "bertini2/flint_version.hpp"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


namespace {

// A ball holding a square root, at `bits` of precision.
struct AcbSqrt
{
    acb_t value;
    AcbSqrt(long radicand, slong bits)
    {
        acb_init(value);
        acb_set_si(value, radicand);
        acb_sqrt(value, value, bits);
    }
    ~AcbSqrt() { acb_clear(value); }
};

} // namespace


BOOST_AUTO_TEST_SUITE(flint)


// A product of balls encloses the exact product, and dividing it back out encloses the factor.
BOOST_AUTO_TEST_CASE(acb_product_and_quotient_enclose_the_exact_values)
{
    const slong working_bits = 64;
    const slong reference_bits = 300;

    AcbSqrt three(3, working_bits), five(5, working_bits);
    AcbSqrt fifteen(15, reference_bits), three_sharp(3, reference_bits);

    acb_t product, quotient;
    acb_init(product); acb_init(quotient);

    acb_mul(product, three.value, five.value, working_bits);
    BOOST_CHECK(acb_contains(product, fifteen.value));
    BOOST_CHECK(!acb_is_exact(product)); // sqrt(3)*sqrt(5) is irrational, so the ball has width

    acb_div(quotient, product, five.value, working_bits);
    BOOST_CHECK(acb_contains(quotient, three_sharp.value));

    acb_clear(product); acb_clear(quotient);
}


// Division by a ball containing zero gives a non-finite ball, not a wrong one.
BOOST_AUTO_TEST_CASE(acb_division_by_a_ball_around_zero_is_not_finite)
{
    const slong bits = 64;
    acb_t numerator, straddling, result;
    acb_init(numerator); acb_init(straddling); acb_init(result);

    acb_set_si(numerator, 7);
    acb_zero(straddling);
    mag_set_ui(arb_radref(acb_realref(straddling)), 2); // [0 +/- 2]

    acb_div(result, numerator, straddling, bits);
    BOOST_CHECK(!acb_is_finite(result));

    acb_clear(numerator); acb_clear(straddling); acb_clear(result);
}


// mpfr -> arf is exact, and so is the way back at the same precision.
BOOST_AUTO_TEST_CASE(arf_round_trips_an_mpfr_value_exactly)
{
    const mpfr_prec_t bits = 150;
    mpfr_t seven_thirds, back;
    mpfr_init2(seven_thirds, bits);
    mpfr_init2(back, bits);
    mpfr_set_ui(seven_thirds, 7, MPFR_RNDN);
    mpfr_div_ui(seven_thirds, seven_thirds, 3, MPFR_RNDN);

    arf_t a;
    arf_init(a);
    arf_set_mpfr(a, seven_thirds);
    const int inexact = arf_get_mpfr(back, a, MPFR_RNDN);

    BOOST_CHECK_EQUAL(inexact, 0);
    BOOST_CHECK(mpfr_equal_p(seven_thirds, back));

    arf_clear(a);
    mpfr_clear(seven_thirds); mpfr_clear(back);
}


// A 337-bit integer crosses from b2's GMP into FLINT and back unchanged: FLINT and b2 agree on
// the mpz layout, which they would not if two GMPs were loaded.  The versions are printed so a
// CI log says which libraries each platform linked.
BOOST_AUTO_TEST_CASE(flint_and_bertini_share_one_gmp_and_one_mpfr)
{
    std::cout << "gmp " << gmp_version << ", mpfr " << mpfr_get_version()
              << " (header " << MPFR_VERSION_STRING << "), flint " << FLINT_VERSION
              << " (runtime " << flint_version << ")\n";
    BOOST_CHECK_EQUAL(std::string(mpfr_get_version()), std::string(MPFR_VERSION_STRING));
    BOOST_CHECK_EQUAL(std::string(flint_version), std::string(FLINT_VERSION));
    BOOST_CHECK_EQUAL(bertini::FlintVersion(), std::string(FLINT_VERSION)); // the library links FLINT too

    mpz_t big, back;
    mpz_init(big); mpz_init(back);
    mpz_ui_pow_ui(big, 7, 120);

    fmpz_t f;
    fmpz_init(f);
    fmpz_set_mpz(f, big);
    fmpz_get_mpz(back, f);
    BOOST_CHECK_EQUAL(mpz_cmp(big, back), 0);

    fmpz_clear(f);
    mpz_clear(big); mpz_clear(back);
}


// numpy zero-fills new buffers of a dtype that needs initialization (ADR-0006).  For MPFR that
// zero is an uninitialized value and reading it crashes; for Arb it is an exact zero, usable and
// clearable as it stands.
BOOST_AUTO_TEST_CASE(a_zero_filled_acb_is_an_exact_zero)
{
    acb_t z;
    std::memset(z, 0, sizeof(acb_struct));
    BOOST_CHECK(acb_is_zero(z));
    BOOST_CHECK(acb_is_exact(z));

    acb_set_si(z, 11);
    acb_mul(z, z, z, 64);
    BOOST_CHECK(acb_equal_si(z, 121));
    acb_clear(z);
}


namespace {

using RealInterval = boost::numeric::interval<double>;
using UnprotectedRealInterval = boost::numeric::interval_lib::unprotect<RealInterval>::type;

// A complex rectangle over a real interval type, by the textbook formulas: every operation is
// the interval library's, none of the rounding is ours.
template <typename IntervalT>
struct Rectangle
{
    IntervalT re, im;
};

template <typename IntervalT>
Rectangle<IntervalT> operator*(Rectangle<IntervalT> const& a, Rectangle<IntervalT> const& b)
{
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

template <typename IntervalT>
Rectangle<IntervalT> operator+(Rectangle<IntervalT> const& a, Rectangle<IntervalT> const& b)
{
    return {a.re + b.re, a.im + b.im};
}

constexpr int polynomial_degree = 10;

// p(z) = sum_k c_k z^k with c_k = (k+2) + (k+3) i, by Horner's rule.
template <typename T>
T horner(T const& z, std::vector<T> const& coefficients)
{
    T acc = coefficients[polynomial_degree];
    for (int k = polynomial_degree - 1; k >= 0; --k)
        acc = acc * z + coefficients[static_cast<std::size_t>(k)];
    return acc;
}

template <typename IntervalT>
std::vector<Rectangle<IntervalT>> rectangle_coefficients()
{
    std::vector<Rectangle<IntervalT>> c;
    for (int k = 0; k <= polynomial_degree; ++k)
        c.push_back({IntervalT(double(k + 2)), IntervalT(double(k + 3))});
    return c;
}

const double z_real = 0.7, z_imag = 0.3;

// The same Horner in acb, through the in-place calls a wrapper would use.
void acb_horner(acb_t acc, acb_t const z, acb_ptr coefficients, slong bits)
{
    acb_set(acc, coefficients + polynomial_degree);
    for (int k = polynomial_degree - 1; k >= 0; --k)
    {
        acb_mul(acc, acc, z, bits);
        acb_add(acc, acc, coefficients + k, bits);
    }
}

} // namespace


// The Boost rectangle, evaluated with the rounding mode set once for the whole evaluation,
// encloses Arb's far tighter enclosure of the same polynomial value, and has width: rounding
// happened in both directions, so the optimizer respected the rounding mode.
BOOST_AUTO_TEST_CASE(boost_interval_rectangle_encloses_the_arb_evaluation)
{
    const slong bits = 256;
    acb_t z, value;
    acb_init(z); acb_init(value);
    acb_set_d_d(z, z_real, z_imag);
    acb_ptr coefficients = _acb_vec_init(polynomial_degree + 1);
    for (int k = 0; k <= polynomial_degree; ++k)
        acb_set_si_si(coefficients + k, k + 2, k + 3);
    acb_horner(value, z, coefficients, bits);

    arf_t lo, hi;
    arf_init(lo); arf_init(hi);
    auto bounds = [&](arb_t const part, double& lower, double& upper) {
        arb_get_lbound_arf(lo, part, bits);
        arb_get_ubound_arf(hi, part, bits);
        lower = arf_get_d(lo, ARF_RND_FLOOR);
        upper = arf_get_d(hi, ARF_RND_CEIL);
    };
    double re_lower, re_upper, im_lower, im_upper;
    bounds(acb_realref(value), re_lower, re_upper);
    bounds(acb_imagref(value), im_lower, im_upper);

    Rectangle<UnprotectedRealInterval> rect;
    {
        RealInterval::traits_type::rounding rounding_scope;
        auto coefficients_rect = rectangle_coefficients<UnprotectedRealInterval>();
        Rectangle<UnprotectedRealInterval> zr{UnprotectedRealInterval(z_real), UnprotectedRealInterval(z_imag)};
        rect = horner(zr, coefficients_rect);
    }

    std::cout << "p(z) real part: arb [" << std::setprecision(17) << re_lower << ", " << re_upper
              << "], boost [" << rect.re.lower() << ", " << rect.re.upper() << "]\n";
    BOOST_CHECK(rect.re.lower() <= re_lower);
    BOOST_CHECK(re_upper <= rect.re.upper());
    BOOST_CHECK(rect.im.lower() <= im_lower);
    BOOST_CHECK(im_upper <= rect.im.upper());
    BOOST_CHECK(rect.re.lower() < rect.re.upper());

    arf_clear(lo); arf_clear(hi);
    _acb_vec_clear(coefficients, polynomial_degree + 1);
    acb_clear(z); acb_clear(value);
}


BOOST_AUTO_TEST_SUITE_END()
