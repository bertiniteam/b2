//This file is part of Bertini 2.
//
//test/blackbox/parsing.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/blackbox/parsing.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/blackbox/parsing.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame

#include <boost/test/unit_test.hpp>

#include "bertini2/system/precon.hpp"
#include "bertini2/blackbox/global_configs.hpp"
#include "bertini2/io/parsing.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"


BOOST_AUTO_TEST_SUITE(blackbox_test)

BOOST_AUTO_TEST_SUITE(parsing_configs)

using namespace bertini;
using mpfr = bertini::mpfr_complex;
using dbl = bertini::dbl;
 
BOOST_AUTO_TEST_CASE(parse1)
{
	
	using AllConfsD = blackbox::config::Configs::All<dbl>::type;
	using AllConfsMP = blackbox::config::Configs::All<mpfr>::type;

std::string config = 
R"(outputlevel: 0;
randomseed: 72;
tracktype: 1;
odepredictor: 8;
finaltol: 1e-12;
endgamenum: 2;
numsamplepoints: 8;
samplefactor: 0.8;
maxcyclenum: 6;
ampmaxprec: 3096;
ampsafetydigits1: 0;
ampsafetydigits2: 0;
maxstepsize: 0.05;
maxnumbersteps: 20000;
securitylevel: 1;
EndpointFiniteThreshold: 1e10;
pathtruncationthreshold: 1e10;
endgamebdry: 0.001;
stepsforincrease: 10;
stepfailfactor: 0.45;
stepsuccessfactor: 1.5;
functiontolerance: 1e-7;
tracktolbeforeeg: 1e-8;
tracktolduringeg: 1e-8;
sharpendigits: 60;
condnumthreshold: 1e300;
maxstepsbeforenewton: 0;
maxnewtonits: 1;)";

	auto results_double = bertini::parsing::classic::ConfigParser<AllConfsD>::Parse(config);
	auto results_mp = bertini::parsing::classic::ConfigParser<AllConfsMP>::Parse(config);

	BOOST_CHECK_EQUAL(std::get<algorithm::RandomConfig>(results_double).random_seed, 72ul);
	BOOST_CHECK_EQUAL(std::get<algorithm::RandomConfig>(results_mp).random_seed, 72ul);
}



BOOST_AUTO_TEST_SUITE_END() // end the parsing sub-suite



BOOST_AUTO_TEST_SUITE(parser_errors)

BOOST_AUTO_TEST_CASE(system_missing_semicolon)
{
	// "variable_group x, y" is missing a semicolon — expectation operator fires
	std::string bad = "variable_group x, y\nfunction f;\nf = x+y;";
	BOOST_CHECK_THROW(bertini::System{bad}, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(system_syntax_error_in_expression)
{
	std::string bad = "variable_group x, y;\nfunction f;\nf = x + * y;";
	BOOST_CHECK_THROW(bertini::System{bad}, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(system_garbage_input)
{
	// Completely nonsensical input — parser cannot make progress
	BOOST_CHECK_THROW(bertini::System{"@#$% not bertini at all"}, std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END() // end parser_errors suite



BOOST_AUTO_TEST_SUITE(unary_minus_precedence)

using dbl = bertini::dbl;

namespace {
	// parse "f = <expr>" over variable_group x,y,z and evaluate at the given point
	dbl ParseEval(std::string const& expr, dbl x, dbl y, dbl z)
	{
		bertini::System s{"variable_group x, y, z;\nfunction f;\nf = " + expr + ";"};
		bertini::Vec<dbl> pt(3);
		pt << x, y, z;
		return s.Eval(pt)(0);
	}
}

// Regression for the grammar bug where a leading unary '-' negated the ENTIRE following
// expression ("-y+x" parsed as "-(y+x)") instead of just its operand.  Unary +/- bind one
// factor_, so "-y+x" is (-y)+x and "-x^2" is -(x^2).
BOOST_AUTO_TEST_CASE(leading_minus_negates_only_its_operand)
{
	const dbl x(2,0), y(5,0), z(3,0);
	// "-y+x" == x-y == -3, NOT -(y+x) == -7
	BOOST_CHECK_SMALL(std::abs(ParseEval("-y+x", x,y,z) - ParseEval("x-y", x,y,z)), 1e-12);
	BOOST_CHECK_SMALL(std::abs(ParseEval("-y+x", x,y,z) - dbl(-3,0)),               1e-12);
}

BOOST_AUTO_TEST_CASE(leading_minus_on_a_parenthesized_sum)
{
	const dbl x(2,0), y(5,0), z(3,0);
	// the bug found via round-tripping: "-(y-z)+x" == x-(y-z) == 0, NOT -((y-z)+x)
	BOOST_CHECK_SMALL(std::abs(ParseEval("-(y-z)+x", x,y,z) - ParseEval("x-(y-z)", x,y,z)), 1e-12);
	BOOST_CHECK_SMALL(std::abs(ParseEval("-(y-z)+x", x,y,z) - dbl(0,0)),                    1e-12);
}

BOOST_AUTO_TEST_CASE(all_negative_sum_is_unchanged)
{
	const dbl x(2,0), y(5,0), z(3,0);
	// "-y-x" == -(y+x) == -7 (here the greedy reading happened to agree)
	BOOST_CHECK_SMALL(std::abs(ParseEval("-y-x", x,y,z) - dbl(-7,0)), 1e-12);
}

BOOST_AUTO_TEST_CASE(unary_minus_binds_looser_than_power)
{
	const dbl x(2,0), y(5,0), z(3,0);
	// "-x^2" == -(x^2) == -4, NOT (-x)^2 == 4
	BOOST_CHECK_SMALL(std::abs(ParseEval("-x^2", x,y,z) - dbl(-4,0)), 1e-12);
}

BOOST_AUTO_TEST_CASE(unary_minus_then_product)
{
	const dbl x(2,0), y(5,0), z(3,0);
	// "-x*y" == -(x*y) == -10
	BOOST_CHECK_SMALL(std::abs(ParseEval("-x*y", x,y,z) - dbl(-10,0)), 1e-12);
}

BOOST_AUTO_TEST_SUITE_END() // unary_minus_precedence



BOOST_AUTO_TEST_SUITE_END() // end the blackbox suite
