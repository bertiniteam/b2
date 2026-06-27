//This file is part of Bertini 2.
//
//system_printing_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_printing_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file system_printing_test.cpp

\brief Human-facing printing of block-composed systems (System::Describe / operator<<): terse
placeholders by default, actual numbers + underlying functions when verbose, and none of the old
debugging noise.
*/

#include <sstream>
#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/block.hpp"

BOOST_AUTO_TEST_SUITE(system_printing_suite)

using namespace bertini;
using bertini::node::Variable;

static std::string Terse(System const& s)   { std::ostringstream ss; s.Describe(ss, false); return ss.str(); }
static std::string Verbose(System const& s) { std::ostringstream ss; s.Describe(ss, true);  return ss.str(); }
static bool Has(std::string const& hay, std::string const& needle) { return hay.find(needle) != std::string::npos; }

// none of the old debugging leftovers should appear.
static void NoNoise(std::string const& t)
{
	BOOST_CHECK(!Has(t, "current variable values"));
	BOOST_CHECK(!Has(t, "differentiated"));
	BOOST_CHECK(!Has(t, "unnamed_function"));
}


BOOST_AUTO_TEST_CASE(plain_polynomial)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System s; s.AddVariableGroup(VariableGroup{x, y});
	s.AddFunction(x*x + y*y - node::Integer::Make(1));
	s.AddFunction(x - y);

	std::string t = Terse(s);
	BOOST_CHECK(Has(t, "f_0 = "));
	BOOST_CHECK(Has(t, "f_1 = "));
	NoNoise(t);
}

BOOST_AUTO_TEST_CASE(randomization_placeholder_and_underlying)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System o; o.AddVariableGroup(VariableGroup{x, y});
	o.AddFunction(x*x + y*y - node::Integer::Make(1));
	o.AddFunction(x*y);
	o.AddFunction(x*x + y*y - x - y);
	System r = o.Randomize();

	std::string t = Terse(r);
	BOOST_CHECK(Has(t, "R . g"));
	BOOST_CHECK(Has(t, "(R: 2x3 randomization matrix)"));
	BOOST_CHECK(Has(t, "g_0 = "));
	BOOST_CHECK(Has(t, "g_2 = "));
	BOOST_CHECK(!Has(t, "R =\n"));               // the matrix is not in the terse form
	NoNoise(t);

	std::string v = Verbose(r);
	BOOST_CHECK(Has(v, "R ="));                  // ... but it is in verbose
}

BOOST_AUTO_TEST_CASE(linear_forms_placeholder_vs_actual)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System m; m.AddVariableGroup(VariableGroup{x, y});
	m.AddFunction(x*x + y*y - node::Integer::Make(1));        // row 0: polynomial
	Mat<mpfr_complex> M(1, 3); M << mpfr_complex(2), mpfr_complex(1), mpfr_complex(-1);
	m.AddBlock(blocks::LinearFormsBlock(2, M));               // row 1: 2x + y - 1

	std::string t = Terse(m);
	BOOST_CHECK(Has(t, "f_0 = "));
	BOOST_CHECK(Has(t, "f_1 = c.[x, y, 1]"));    // both rows visible; structured row is a placeholder

	std::string v = Verbose(m);
	BOOST_CHECK(Has(v, "*x") && Has(v, "*y"));   // actual coefficients shown
}

BOOST_AUTO_TEST_CASE(moving_homotopy_blend)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y}); fixed.AddFunction(x*x + y*y - node::Integer::Make(1));
	System sm; sm.AddVariableGroup(VariableGroup{x, y}); sm.AddFunction(y);
	System em; em.AddVariableGroup(VariableGroup{x, y}); em.AddFunction(y - x);
	System H = MakeMovingHomotopy(fixed, sm, em, "t", node::Complex::Make(mpfr_complex("0.6", "0.8")));

	std::string t = Terse(H);
	BOOST_CHECK(Has(t, "f_0 = "));               // the fixed polynomial row
	BOOST_CHECK(Has(t, "blend of 2 systems"));
	BOOST_CHECK(Has(t, "path variable: t"));
	BOOST_CHECK(!Has(t, "f_1..f_1"));            // single moving row reads f_1

	std::string v = Verbose(H);
	BOOST_CHECK(Has(v, "A_0 = ") && Has(v, "B_0 = "));   // operand functions listed
}

BOOST_AUTO_TEST_SUITE_END()
