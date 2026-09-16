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
using real_mp  = bertini::real_mp;
using complex_mp = bertini::complex_mp;
using mpq_rational = bertini::mpq_rational;

// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------

// Simulate a stale non-local-static real_mp: created at STARTUP_PREC
// then used in arithmetic at a different target precision.
static constexpr unsigned STARTUP_PREC = 20;  // BMP default at program start

// -----------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(precision_propagation)
// -----------------------------------------------------------------------


// --- 1. DefaultPrecision is uniform across float and complex ------------

BOOST_AUTO_TEST_CASE(default_precision_sets_both_float_and_complex)
{
	DefaultPrecision(30);
	BOOST_CHECK_EQUAL(real_mp::default_precision(),  30u);
	BOOST_CHECK_EQUAL(complex_mp::default_precision(), 30u);

	DefaultPrecision(16);
	BOOST_CHECK_EQUAL(real_mp::default_precision(),  16u);
	BOOST_CHECK_EQUAL(complex_mp::default_precision(), 16u);

	DefaultPrecision(50);
	BOOST_CHECK_EQUAL(real_mp::default_precision(),  50u);
	BOOST_CHECK_EQUAL(complex_mp::default_precision(), 50u);
}


// --- 2. Precision policy is uniform after DefaultPrecision() -----------
//  (only meaningful with expression templates on, but harmless either way)

#ifdef BMP_EXPRESSION_TEMPLATES
BOOST_AUTO_TEST_CASE(default_precision_sets_preserve_related_on_both_types)
{
	using bmp = boost::multiprecision::variable_precision_options;
	DefaultPrecision(30);
	BOOST_CHECK(real_mp::thread_default_variable_precision_options()
	            == bmp::preserve_related_precision);
	BOOST_CHECK(complex_mp::thread_default_variable_precision_options()
	            == bmp::preserve_related_precision);
}
#endif


// --- 3. mpq_rational carries no MPFR precision state ------------------
//
// Key property: mpq_rational * real_mp(prec=N) → result at prec=N.
// This is the fix for EndgameConfig::sample_factor.

BOOST_AUTO_TEST_CASE(mpq_rational_times_mpfr_float_result_at_target_precision)
{
	DefaultPrecision(16);
	real_mp t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	mpq_rational half{1, 2};
	real_mp result = t * half;
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(mpq_rational_times_mpfr_float_at_high_precision)
{
	DefaultPrecision(50);
	real_mp t("0.3");
	mpq_rational r{647, 1000};
	real_mp result = t * r;
	BOOST_CHECK_EQUAL(result.precision(), 50u);
}


// --- 4. Stale real_mp at prec > target contaminates arithmetic -----
//
// DOCUMENTARY test: demonstrates the PROBLEM that motivates using
// mpq_rational for config fields.  With preserve_related_precision,
// an real_mp stale at prec=STARTUP_PREC poisons results when
// STARTUP_PREC > target precision.

BOOST_AUTO_TEST_CASE(stale_mpfr_float_contaminates_when_startup_prec_exceeds_target)
{
	// Create a stale config value at "startup precision" (20 digits).
	DefaultPrecision(STARTUP_PREC);
	real_mp stale_factor("0.5");
	BOOST_REQUIRE_EQUAL(stale_factor.precision(), STARTUP_PREC);

	// Now switch to a lower target precision.
	DefaultPrecision(16);
	real_mp t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	// With preserve_related_precision the result takes max(stale, target) = 20 — NOT 16.
	// This test DOCUMENTS the contamination and is expected to fail/warn until
	// all config fields that appear in non-local statics are changed to mpq_rational.
#ifdef BMP_EXPRESSION_TEMPLATES
	real_mp contaminated = t * stale_factor;
	BOOST_WARN_EQUAL(contaminated.precision(), 16u);  // WARN: currently 20 (contaminated)
	BOOST_CHECK_EQUAL(contaminated.precision(), STARTUP_PREC);  // confirms contamination
#endif
}


// --- 5. EndgameConfig::sample_factor is mpq_rational (regression guard)
//
// If someone changes sample_factor back to real_mp, this test
// catches the precision contamination at prec < STARTUP_PREC.

BOOST_AUTO_TEST_CASE(endgame_config_sample_factor_no_precision_contamination)
{
	// Use the static DefaultConstruct value — exactly what happens in production
	// when an endgame is constructed without an explicit config.
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::endgame::EndgameConfig>::value;

	DefaultPrecision(16);
	real_mp t("0.3");
	BOOST_REQUIRE_EQUAL(t.precision(), 16u);

	real_mp result = t * real_mp(cfg.sample_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(endgame_config_default_sample_factor_is_one_half)
{
	bertini::endgame::EndgameConfig cfg;
	BOOST_CHECK_EQUAL(cfg.sample_factor, mpq_rational(1, 2));
}


// --- 6. SteppingConfig step-size fields should not contaminate ---------
//
// These tests currently FAIL because SteppingConfig uses real_mp for
// its step-size fields — they get initialized at prec=STARTUP_PREC in
// DefaultConstruct<SteppingConfig>::value.  They will PASS once those
// fields are changed to mpq_rational (or double for the min_step_size).

BOOST_AUTO_TEST_CASE(stepping_config_success_factor_no_precision_contamination)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	real_mp step("0.1");
	BOOST_REQUIRE_EQUAL(step.precision(), 16u);

	real_mp result = step * real_mp(cfg.step_size_success_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(stepping_config_fail_factor_no_precision_contamination)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	real_mp step("0.1");
	BOOST_REQUIRE_EQUAL(step.precision(), 16u);

	real_mp result = step * real_mp(cfg.step_size_fail_factor);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}

BOOST_AUTO_TEST_CASE(stepping_config_max_step_size_no_stale_precision)
{
	auto const& cfg = bertini::detail::DefaultConstruct<bertini::tracking::SteppingConfig>::value;

	DefaultPrecision(16);
	// max_step_size itself should be at (or convertible to) target precision without contamination
	real_mp result = real_mp(cfg.max_step_size);
	BOOST_CHECK_EQUAL(result.precision(), 16u);
}


// -----------------------------------------------------------------------


// --- SLP precision follows its EVALUATION ARGUMENTS (#377) -------------
//
// The compiled Program is a precision-independent tape of operations; only the Memory
// holding values carries digits.  So the Memory's precision is an artifact of the current
// evaluation, never an invariant, and evaluation must re-tag rather than refuse.
//
// This used to throw, and -- worse -- there was no way out of it.  A Memory takes its
// precision from the ambient DefaultPrecision() when the program is lazily compiled, while
// the owning System keeps whatever it was told, so the two diverge the moment anything
// moves the ambient default (an AMP tracker or endgame does, routinely).  Both
// System::precision() and StraightLineProgram::precision() short-circuit when handed the
// value they already hold, so neither could repair the divergence: the System was wedged.

namespace {
	/// f(x,y) = x^2 - y, whose Jacobian [2x, -1] is exact at any precision.
	bertini::System TwoVarSystem()
	{
		using bertini::node::Variable;
		auto x = Variable::Make("x"), y = Variable::Make("y");
		bertini::System sys;
		bertini::VariableGroup vars{x, y};
		sys.AddVariableGroup(vars);
		sys.AddFunction(pow(x, 2) - y);
		return sys;
	}

	bertini::Vec<complex_mp> PointAt(unsigned prec, std::string const& re)
	{
		auto saved = DefaultPrecision();
		DefaultPrecision(prec);
		bertini::Vec<complex_mp> p(2);
		p << complex_mp(real_mp(re), real_mp("0")), complex_mp(real_mp("1"), real_mp("0"));
		Precision(p, prec);
		DefaultPrecision(saved);
		return p;
	}
}

BOOST_AUTO_TEST_CASE(slp_evaluates_after_the_ambient_default_moved_underneath_it)
{
	DefaultPrecision(30);
	auto sys = TwoVarSystem();
	sys.precision(30);

	// an AMP tracker or endgame leaves the ambient default somewhere else entirely
	DefaultPrecision(16);

	// evaluating at the system's OWN precision must work.  before the fix this threw
	// "variable_values and SLP must be of same precision.  respective precisions: 30 16"
	auto pt = PointAt(30, "3");
	BOOST_CHECK_NO_THROW(sys.Eval(pt));
	auto v = sys.Eval(pt);
	BOOST_CHECK_EQUAL(v.size(), 1);
	BOOST_CHECK(abs(v(0) - complex_mp(8)) < 1e-25);   // 3^2 - 1 = 8

	BOOST_CHECK_NO_THROW(sys.Jacobian(pt));
	auto J = sys.Jacobian(pt);
	BOOST_CHECK(abs(J(0,0) - complex_mp(6)) < 1e-25); // d/dx = 2x = 6
	BOOST_CHECK(abs(J(0,1) + complex_mp(1)) < 1e-25); // d/dy = -1
}

BOOST_AUTO_TEST_CASE(slp_retags_when_evaluated_at_a_sequence_of_precisions)
{
	DefaultPrecision(16);
	auto sys = TwoVarSystem();

	// the same system, evaluated up and down a ladder of precisions, in one lifetime
	for (unsigned prec : {16u, 50u, 30u, 100u, 20u})
	{
		auto pt = PointAt(prec, "3");
		BOOST_CHECK_NO_THROW(sys.Eval(pt));
		auto v = sys.Eval(pt);
		BOOST_CHECK(abs(v(0) - complex_mp(8)) < 1e-15);
	}
}

BOOST_AUTO_TEST_CASE(slp_retag_rebuilds_constants_it_does_not_zero_pad_them)
{
	// The mean version.  A constant that is INEXACT in decimal -- 1/3 -- must come back at
	// full accuracy after a re-tag upward, which is only true if re-tagging refills the
	// constants from their exact recipes.  Padding a 16-digit 1/3 out to 100 digits would
	// leave the tail zero and the residual ~1e-17, not ~1e-100.
	using bertini::node::Variable;
	using bertini::node::Rational;
	auto x = Variable::Make("x");
	bertini::System sys;
	bertini::VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddFunction(x - Rational::Make(mpq_rational(1, 3)));

	DefaultPrecision(16);
	{	// force the lazy compile to happen at 16 digits
		auto saved = DefaultPrecision();
		bertini::Vec<complex_mp> p(1);
		p << complex_mp(real_mp("1"), real_mp("0"));
		Precision(p, 16);
		BOOST_CHECK_NO_THROW(sys.Eval(p));
		DefaultPrecision(saved);
	}

	// now evaluate AT 1/3 to 100 digits: the residual must be ~1e-100, not ~1e-17
	DefaultPrecision(100);
	bertini::Vec<complex_mp> third(1);
	third << complex_mp(real_mp(mpq_rational(1, 3)), real_mp("0"));
	Precision(third, 100);
	auto v = sys.Eval(third);
	BOOST_CHECK_MESSAGE(abs(v(0)) < 1e-90,
		"constant was padded rather than rebuilt: residual " << abs(v(0)));
}

BOOST_AUTO_TEST_CASE(slp_memory_ends_at_the_max_of_its_arguments_precisions)
{
	// The path variable arrives AFTER the variables are already in memory, so re-tagging on
	// it may only ever raise -- lowering would truncate the variables just written.
	using bertini::node::Variable;
	auto x = Variable::Make("x");
	auto t = Variable::Make("t");
	bertini::System sys;
	bertini::VariableGroup vars{x};
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(x * t);

	DefaultPrecision(30);
	sys.precision(30);
	auto pt = PointAt(30, "2");
	bertini::Vec<complex_mp> one(1);
	one << pt(0);                      // x = 2 at 30 digits

	DefaultPrecision(60);
	complex_mp time(real_mp("3"), real_mp("0"));
	Precision(time, 60);

	BOOST_CHECK_NO_THROW(sys.Eval(one, time));
	auto v = sys.Eval(one, time);
	BOOST_CHECK(abs(v(0) - complex_mp(6)) < 1e-25);
	BOOST_CHECK_GE(Precision(v(0)), 30u);
}

BOOST_AUTO_TEST_SUITE_END()
