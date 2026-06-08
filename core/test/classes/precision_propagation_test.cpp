//This file is part of Bertini 2.
//
//test/classes/precision_propagation_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/classes/precision_propagation_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/classes/precision_propagation_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

#include <boost/test/unit_test.hpp>

#include "bertini2/bertini.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/detail/enable_permuted_arguments.hpp"

using bertini::DefaultPrecision;
using bertini::Precision;
using mpfr_float  = bertini::mpfr_float;
using mpfr_complex = bertini::mpfr_complex;
using mpq_rational = bertini::mpq_rational;

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

// Simulate a stale non-local-static mpfr_float: created at STARTUP_PREC
// then used in arithmetic at a different target precision.
static constexpr unsigned STARTUP_PREC = 20;  // BMP default at program start

// -----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(precision_propagation)
// -----------------------------------------------------------------------


// --- 1. DefaultPrecision is uniform across float and complex ------------

BOOST_AUTO_TEST_CASE(default_precision_sets_both_float_and_complex)
{
	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(mpfr_float::default_precision(),  30u);
	BOOST_CHECK_EQUAL(mpfr_complex::default_precision(), 30u);

	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(mpfr_float::default_precision(),  16u);
	BOOST_CHECK_EQUAL(mpfr_complex::default_precision(), 16u);

	DefaultPrecision(50);
	BOOST_CHECK_EQUAL(mpfr_float::default_precision(),  50u);
	BOOST_CHECK_EQUAL(mpfr_complex::default_precision(), 50u);
}


// --- 2. Precision policy is uniform after DefaultPrecision() -----------
//  (only meaningful with expression templates on, but harmless either way)

#ifdef BMP_EXPRESSION_TEMPLATES
BOOST_AUTO_TEST_CASE(default_precision_sets_preserve_related_on_both_types)
{
	using bmp = boost::multiprecision::variable_precision_options;
	DefaultPrecision(30);
	BOOST_CHECK(mpfr_float::thread_default_variable_precision_options()
	            == bmp::preserve_related_precision);
	BOOST_CHECK(mpfr_complex::thread_default_variable_precision_options()
	            == bmp::preserve_related_precision);
}
#endif


// --- 3. mpq_rational carries no MPFR precision state ------------------
//
// Key property: mpq_rational * mpfr_float(prec=N) → result at prec=N.
// This is the fix for EndgameConfig::sample_factor.

BOOST_AUTO_TEST_CASE(mpq_rational_times_mpfr_float_result_at_target_precision)
{
	DefaultPrecision(16);
	mpfr_float t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	mpq_rational half{1, 2};
	mpfr_float result = t * half;
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(mpq_rational_times_mpfr_float_at_high_precision)
{
	DefaultPrecision(50);
	mpfr_float t("0.3");
	mpq_rational r{647, 1000};
	mpfr_float result = t * r;
	BOOST_CHECK_EQUAL(result.precision(), 50u);
}


// --- 4. Stale mpfr_float at prec > target contaminates arithmetic -----
//
// DOCUMENTARY test: demonstrates the PROBLEM that motivates using
// mpq_rational for config fields.  With preserve_related_precision,
// an mpfr_float stale at prec=STARTUP_PREC poisons results when
// STARTUP_PREC > target precision.

BOOST_AUTO_TEST_CASE(stale_mpfr_float_contaminates_when_startup_prec_exceeds_target)
{
	// Create a stale config value at "startup precision" (20 digits).
	DefaultPrecision(STARTUP_PREC);
	mpfr_float stale_factor("0.5");
	BOOST_REQUIRE_EQUAL(stale_factor.precision(), STARTUP_PREC);

	// Now switch to a lower target precision.
	DefaultPrecision(16);
	mpfr_float t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	// With preserve_related_precision the result takes max(stale, target) = 20 — NOT 16.
	// This test DOCUMENTS the contamination and is expected to fail/warn until
	// all config fields that appear in non-local statics are changed to mpq_rational.
#ifdef BMP_EXPRESSION_TEMPLATES
	mpfr_float contaminated = t * stale_factor;
	BOOST_WARN_EQUAL(contaminated.precision(), 16u);  // WARN: currently 20 (contaminated)
	BOOST_CHECK_EQUAL(contaminated.precision(), STARTUP_PREC);  // confirms contamination
#endif
}


// --- 5. EndgameConfig::sample_factor is mpq_rational (regression guard)
//
// If someone changes sample_factor back to mpfr_float, this test
// catches the precision contamination at prec < STARTUP_PREC.

BOOST_AUTO_TEST_CASE(endgame_config_sample_factor_no_precision_contamination)
{
	// Use the static DefaultConstruct value — exactly what happens in production
	// when an endgame is constructed without an explicit config.
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::endgame::EndgameConfig>::value;

	DefaultPrecision(16);
	mpfr_float t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	mpfr_float result = t * mpfr_float(cfg.sample_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(endgame_config_default_sample_factor_is_one_half)
{
	bertini::endgame::EndgameConfig cfg;
	BOOST_CHECK_EQUAL(cfg.sample_factor, mpq_rational(1, 2));
}


// --- 6. SteppingConfig step-size fields should not contaminate ---------
//
// These tests currently FAIL because SteppingConfig uses mpfr_float for
// its step-size fields — they get initialized at prec=STARTUP_PREC in
// DefaultConstruct<SteppingConfig>::value.  They will PASS once those
// fields are changed to mpq_rational (or double for the min_step_size).

BOOST_AUTO_TEST_CASE(stepping_config_success_factor_no_precision_contamination)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	mpfr_float step("0.1");
	BOOST_REQUIRE_EQUAL(step.precision(), 16u);

	mpfr_float result = step * mpfr_float(cfg.step_size_success_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(stepping_config_fail_factor_no_precision_contamination)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	mpfr_float step("0.1");
	BOOST_REQUIRE_EQUAL(step.precision(), 16u);

	mpfr_float result = step * mpfr_float(cfg.step_size_fail_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(stepping_config_max_step_size_no_stale_precision)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	// max_step_size itself should be at (or convertible to) target precision without contamination
	mpfr_float result = mpfr_float(cfg.max_step_size);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}


// -----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE_END()
