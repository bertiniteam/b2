//This file is part of Bertini 2.
//
//system_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

//  system_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.

/**
\file system_test.cpp Unit testing for the bertini::System class.
*/

#include <boost/test/unit_test.hpp>



#include <sstream>

#include "bertini2/system/system.hpp"
#include "bertini2/system/precon.hpp"
#include "bertini2/system/slice.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

#include "externs.hpp"

namespace utf = boost::unit_test;

using Variable = bertini::node::Variable;


BOOST_AUTO_TEST_SUITE(system_class)

using Var = std::shared_ptr<bertini::node::Variable>;

using mpfr = bertini::complex_mp;

using namespace bertini;
/**
\class bertini::System
\test \b system_make_a_system_at_all Confirms that can default construct a System.
*/
BOOST_AUTO_TEST_CASE(system_make_a_system_at_all)
{
	System S;
}


/**
\class bertini::System
\test \b system_create_parser Confirms that can parse a system defined by a Bertini Classic style string.
*/
BOOST_AUTO_TEST_CASE(system_create_parser)
{
	System sys;
	std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";

	bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	BOOST_CHECK(s);
	BOOST_CHECK(!sys.IsHomogeneous());

}


/**
\class bertini::System
\test \b parsed_system_has_functions_and_evaluates A parsed system must actually contain its
functions and evaluate them.  (Regression: when classic functions became eager-bound, a parse
path that skipped the post-parse emit produced a silently empty, size-0 system -- the earlier
parse tests only checked parse success / homogeneity, never the function count or a value.)
*/
BOOST_AUTO_TEST_CASE(parsed_system_has_functions_and_evaluates)
{
	System sys;
	std::string str = "variable_group x, y; function f, g; f = x*y; g = x + y;";
	bool ok = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	BOOST_CHECK(ok);
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2u);

	Vec<complex_dbl> pt(2); pt << complex_dbl(2,0), complex_dbl(3,0);
	auto v = sys.Eval(pt);
	BOOST_REQUIRE_EQUAL(v.size(), 2);
	BOOST_CHECK_SMALL(std::abs(v(0) - complex_dbl(6,0)), 1e-12);   // x*y at (2,3)
	BOOST_CHECK_SMALL(std::abs(v(1) - complex_dbl(5,0)), 1e-12);   // x+y at (2,3)
}


/**
\class bertini::System
\test \b parsed_subfunction_is_a_named_expression A `s = expr` subfunction parses to an immutable
NamedExpression embedded in the function that references it; the System evaluates correctly, prints
the function with `s` by name, and lists the named subexpression's value below (discovered, not stored).
*/
BOOST_AUTO_TEST_CASE(parsed_subfunction_is_a_named_expression)
{
	System sys;
	std::string str = "variable_group x, y; function f; s = x*y; f = s + x;";
	bool ok = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	BOOST_CHECK(ok);
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 1u);

	// f = s + x = x*y + x; at (2,3) -> 6 + 2 = 8
	Vec<complex_dbl> pt(2); pt << complex_dbl(2,0), complex_dbl(3,0);
	auto v = sys.Eval(pt);
	BOOST_REQUIRE_EQUAL(v.size(), 1);
	BOOST_CHECK_SMALL(std::abs(v(0) - complex_dbl(8,0)), 1e-12);

	std::stringstream ss; ss << sys;
	const std::string text = ss.str();
	BOOST_CHECK(text.find("named subexpression") != std::string::npos);  // discovered + shown
	BOOST_CHECK(text.find("s = x*y") != std::string::npos);              // its value
}





/**
\class bertini::System
\test \b system_parse_xyz_f1f2_t_pq Confirms that can parse a system defined by a Bertini Classic style string, containing two functions, powers, parameters, a path variable, and three variables.
*/
BOOST_AUTO_TEST_CASE(system_parse_xyz_f1f2_t_pq)
{

	std::string str = "variable_group x, y, z;\n function f1, f2;\n pathvariable t;\n parameter p, q;\n p = t;\n q = 1-t;\n f1 = x*y*z;\n\nf2 = p*q*x - 2^(-5);\n";




	bertini::System sys;
	bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	BOOST_CHECK(s);

}


/**
\class bertini::System
\test \b system_parse_with_subfunctions Confirms that can parse a system with a subfunction.
*/
BOOST_AUTO_TEST_CASE(system_parse_with_subfunctions)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	BOOST_CHECK(s);

	BOOST_CHECK(!sys.IsHomogeneous());
}


/**
\class bertini::System
\test \b system_parse_around_the_unit_circle Tests parsing a system with cosine, sine, \f$\pi\f$, \f$i\f$, and subfunctions.
*/
BOOST_AUTO_TEST_CASE(system_parse_around_the_unit_circle)
{
	std::string str =
 "variable z;\nfunction H;\nparameter q1,q2;\npathvariable t;\nq1 = cos(2*Pi*(1-t));\nq2 = sin(2*Pi*(1-t));\ns = q1 + I*q2;\nH = z^2 - s;\n";



	bertini::System sys;
	bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	BOOST_CHECK(s);

	BOOST_CHECK(!sys.IsHomogeneous());
}



/**
\class bertini::System
\test \b system_parse_around_the_unit_circle_alt Tests parsing an exponential function with time in it as a parameter.
*/
BOOST_AUTO_TEST_CASE(system_parse_around_the_unit_circle_alt)
{
	std::string str = " variable z; function H; parameter s; pathvariable t; s = exp(2*Pi*I*(1-t)); H = z^2 - s; ";


	bertini::System sys;
	bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	BOOST_CHECK(s);

}



//TODO: uncomment this test once error handling has been done for the parsers.
// BOOST_AUTO_TEST_CASE(system_parse_x_y_not_xy)
// {
// 	std::string str = " variable x, y; function f; f = xy;";


// 	bertini::System sys;
// 	std::string::const_iterator iter = str.begin();
// 	std::string::const_iterator end = str.end();
// 	bertini::SystemParser<std::string::const_iterator> S;
// 	bool s = phrase_parse(iter, end, S, boost::spirit::ascii::space, sys);
// 	BOOST_CHECK(!s && iter!=end);

// }



/**
\class bertini::System
\test \b system_differentiate_x Create a system with two functions, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x)
{
	Var x = Variable::Make("x");
	auto f1 = pow(x,2);
	auto f2 = x-1;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<complex_dbl> v(1);
	v << 1.0;

	auto J = S.Jacobian(v);

	BOOST_CHECK_EQUAL(J(0),2.0);
	BOOST_CHECK_EQUAL(J(1),1.0);
}


/**
\class bertini::System
\test \b system_differentiate_x_and_y Create a system with two functions and two variables, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x_and_y)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	auto f1 = pow(x,2)*y/2;
	auto f2 = x-y;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddUngroupedVariable(y);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<complex_dbl> v(2);
	v << 1.0 , 2.0;

	auto J = S.Jacobian(v);

	BOOST_CHECK_THROW(S.Jacobian(v,complex_dbl(0.5)), std::runtime_error);
}


/**
\class bertini::System
\test \b system_differentiate_x_and_t Create a system with two functions, one variable, one time variable, and check that can compute and evaluate its Jacobian.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_x_and_t)
{
	Var x = Variable::Make("x");
	Var t = Variable::Make("t");
	auto f1 = (1-t)*x + t*(1-x);
	auto f2 = x-t;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<complex_dbl> v(1);
	v << 1.0;
	complex_dbl time(0.5,0.2);
	bertini::Mat<complex_dbl> J = S.Jacobian(v,time);

	BOOST_CHECK_THROW(S.Jacobian(v), std::runtime_error);
}




/**
\class bertini::System
\test \b system_homogenize_multiple_variable_groups Homogenize a system with multiple variable groups.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_multiple_variable_groups)
{
	Var x1 = Variable::Make("x1");
	Var x2 = Variable::Make("x2");

	Var y1 = Variable::Make("y1");
	Var y2 = Variable::Make("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*x1 + y1*y2 + x1 + y2 - 1;


	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1);
	S.AddFunction(f2);

	BOOST_CHECK(!S.IsHomogeneous());

	BOOST_CHECK_EQUAL(S.NumHomVariables(),0);

	S.Homogenize();

	BOOST_CHECK(S.IsHomogeneous());
}



/**
\class bertini::System
\test \b system_reorder_by_degree_decreasing For a system with four functions, re-order the functions so they are in decreasing degree.
*/
BOOST_AUTO_TEST_CASE(system_reorder_by_degree_decreasing)
{
	Var x1 = Variable::Make("x1");
	Var x2 = Variable::Make("x2");

	Var y1 = Variable::Make("y1");
	Var y2 = Variable::Make("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*pow(x1,2) + y1*y2 + x1 + y2 - 1;

	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1); // deg 2
	S.AddFunction(f2); // deg 3 
	S.AddFunction(pow(x1,4) + pow(y2,5)); // deg 5
	S.AddFunction(x1 + x2 + y1 + y2); // deg 1


	BOOST_CHECK(!S.IsHomogeneous());
	S.ReorderFunctionsByDegreeDecreasing();

	auto degs = S.Degrees();

	for (auto d = degs.begin(); d != degs.end()-1; d++)
	{
		BOOST_CHECK(*d >= *(d+1));
	}

}



/**
\class bertini::System
\test \b system_reorder_by_degree_increasing For a system with four functions, re-order the functions so they are in increasing degree.
*/
BOOST_AUTO_TEST_CASE(system_reorder_by_degree_increasing)
{
	Var x1 = Variable::Make("x1");
	Var x2 = Variable::Make("x2");

	Var y1 = Variable::Make("y1");
	Var y2 = Variable::Make("y2");


	bertini::VariableGroup v1{x1, x2};
	bertini::VariableGroup v2{y1, y2};

	auto f1 = x1*y1 + x1;
	auto f2 = x2*pow(x1,2) + y1*y2 + x1 + y2 - 1;


	bertini::System S;
	S.AddVariableGroup(v1);
	S.AddVariableGroup(v2);
	
	S.AddFunction(f1);
	S.AddFunction(f2);

	BOOST_CHECK(!S.IsHomogeneous());
	S.ReorderFunctionsByDegreeIncreasing();

	auto degs = S.Degrees();

	for (auto d = degs.begin(); d != degs.end()-1; d++)
	{
		BOOST_CHECK(*d <= *(d+1));
	}

}





/**
\class bertini::System
\test \b system_evaluate_double Evaluate a system in double precision.
*/
BOOST_AUTO_TEST_CASE(system_evaluate_double)
{

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	[[maybe_unused]] bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Vec<complex_dbl> values(2);

	values(0) = complex_dbl(2.0);
	values(1) = complex_dbl(3.0);

	Vec<complex_dbl> v = sys.Eval(values);

	BOOST_CHECK_EQUAL(v(0), 36.0);


	auto J = sys.Jacobian(values);

	double x1 = 2;
	double x2 = 3;

	BOOST_CHECK_EQUAL(J(0,0), 2*x1*x2*x2);
	BOOST_CHECK_EQUAL(J(0,1), x1*x1*2*x2);
}


/**
\class bertini::System
\test \b system_evaluate_mpfr Evaluate a system in multiple precision.
*/
BOOST_AUTO_TEST_CASE(system_evaluate_mpfr)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";

	bertini::System sys;
	[[maybe_unused]] bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	Vec<mpfr> values(2);

	values(0) = mpfr(2);
	values(1) = mpfr(3);

	Vec<mpfr> v = sys.Eval(values);

	BOOST_CHECK_EQUAL(v(0), mpfr(36));


	auto J = sys.Jacobian(values);

	mpfr x1(2);
	mpfr x2(3);
	
	BOOST_CHECK_EQUAL(J(0,0), mpfr(2.0)*x1*x2*x2);
	BOOST_CHECK_EQUAL(J(0,1), x1*x1*mpfr(2.0)*x2);
}


BOOST_AUTO_TEST_CASE(system_jacobian)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");

	// Modest magnitudes keep the high-degree SLP-vs-analytic rounding comfortably inside the
	// 1e-15 tolerance below.
	complex_dbl a(0.5,  0.25);
	complex_dbl b(0.4, -0.30);
	complex_dbl c(0.6,  0.20);

	System sys;

	sys.AddVariableGroup(VariableGroup{x, y, z});

	sys.AddFunction(pow(x,2)*pow(y,3)*pow(z,4) + 1);
	sys.AddFunction(pow(x,3)*pow(y,4)*pow(z,5) + 4);

	Vec<complex_dbl> v(3);
	v << a, b, c;

	auto J = sys.Jacobian(v);


	BOOST_CHECK_SMALL(J(0,0).real() - (2.*a*pow(b,3)*pow(c,4)).real(),        1e-15);
	BOOST_CHECK_SMALL(J(0,1).real() - (3.*pow(a,2)*pow(b,2)*pow(c,4)).real(), 1e-15);
	BOOST_CHECK_SMALL(J(0,2).real() - (4.*pow(a,2)*pow(b,3)*pow(c,3)).real(), 1e-15);

	BOOST_CHECK_SMALL(J(1,0).real() - (3.*pow(a,2)*pow(b,4)*pow(c,5)).real(), 1e-15);
	BOOST_CHECK_SMALL(J(1,1).real() - (4.*pow(a,3)*pow(b,3)*pow(c,5)).real(), 1e-15);
	BOOST_CHECK_SMALL(J(1,2).real() - (5.*pow(a,3)*pow(b,4)*pow(c,4)).real(), 1e-15);

	BOOST_CHECK_SMALL(J(0,0).imag() - (2.*a*pow(b,3)*pow(c,4)).imag(),        1e-15);
	BOOST_CHECK_SMALL(J(0,1).imag() - (3.*pow(a,2)*pow(b,2)*pow(c,4)).imag(), 1e-15);
	BOOST_CHECK_SMALL(J(0,2).imag() - (4.*pow(a,2)*pow(b,3)*pow(c,3)).imag(), 1e-15);

	BOOST_CHECK_SMALL(J(1,0).imag() - (3.*pow(a,2)*pow(b,4)*pow(c,5)).imag(), 1e-15);
	BOOST_CHECK_SMALL(J(1,1).imag() - (4.*pow(a,3)*pow(b,3)*pow(c,5)).imag(), 1e-15);
	BOOST_CHECK_SMALL(J(1,2).imag() - (5.*pow(a,3)*pow(b,4)*pow(c,4)).imag(), 1e-15);


}


/**
\class bertini::System
\test \b add_two_systems Test the arithmetic sum of two Systems.
*/
BOOST_AUTO_TEST_CASE(add_two_systems)
{
	bertini::System sys1, sys2;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys1.AddVariableGroup(vars);
	sys1.AddFunction(y+1);
	sys1.AddFunction(x*y);

	sys2.AddVariableGroup(vars);
	sys2.AddFunction(-y-1);
	sys2.AddFunction(-x*y);

	sys1+=sys2;


	Vec<complex_dbl> values(2);

	values << complex_dbl(2.0), complex_dbl(3.0);

	auto v = sys1.Eval(values);

	BOOST_CHECK_EQUAL(v(0), 0.0);
	BOOST_CHECK_EQUAL(v(1), 0.0);

	auto deg = sys1.Degrees();

	BOOST_CHECK_EQUAL(deg.size(),2);
	if (deg.size()==2)
	{
		BOOST_CHECK_EQUAL(deg[0],1);
		BOOST_CHECK_EQUAL(deg[1],2);
	}


}


/**
\class bertini::System
\test \b add_two_systems_evaluated_in_mpfr Sibling of add_two_systems, but
exercises the Vec<mpfr> evaluation path. This is the C++ analog of the
Python test_add_systems that surfaced the macos-15-intel SIGABRT — the C++
suite previously only covered the Vec<complex_dbl> path.
*/
BOOST_AUTO_TEST_CASE(add_two_systems_evaluated_in_mpfr)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::System sys1, sys2;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys1.AddVariableGroup(vars);
	sys1.AddFunction(y+1);
	sys1.AddFunction(x*y);

	sys2.AddVariableGroup(vars);
	sys2.AddFunction(-y-1);
	sys2.AddFunction(-x*y);

	sys1+=sys2;

	Vec<mpfr> values(2);
	values << mpfr(2), mpfr(3);

	auto v = sys1.Eval(values);

	BOOST_CHECK_EQUAL(v(0), mpfr(0));
	BOOST_CHECK_EQUAL(v(1), mpfr(0));
}


/**
\class bertini::System
\test \b add_system_to_self_doubles_under_function_tree_eval Self-add of a
System should double every function. Forces FunctionTree eval to bypass any
SLP caching, isolating the question "does the symbolic mutation in
operator+= work for an aliased rhs?" to just the function-tree level.
*/
BOOST_AUTO_TEST_CASE(add_system_to_self_doubles_under_function_tree_eval)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys.AddVariableGroup(vars);
	sys.AddFunction(y+1);
	sys.AddFunction(x*y);

	Vec<complex_dbl> values(2);
	values << complex_dbl(2.0), complex_dbl(3.0);
	auto before = sys.Eval(values);   // [4, 6]

	sys += sys;                       // self-add

	auto after = sys.Eval(values);    // expected [8, 12]
	BOOST_CHECK_EQUAL(after(0), 2.0 * before(0));
	BOOST_CHECK_EQUAL(after(1), 2.0 * before(1));
}


/**
\class bertini::System
\test \b operator_plus_equals_invalidates_slp_cache System::operator+=
mutates each Function's entry_node_ via SetRoot. The cached
StraightLineProgram (the default eval method) is built lazily on first
Differentiate(); if we don't invalidate is_differentiated_ here, a
subsequent Eval reads stale values. Pin this so it doesn't regress.
*/
BOOST_AUTO_TEST_CASE(operator_plus_equals_invalidates_slp_cache)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys.AddVariableGroup(vars);
	sys.AddFunction(y+1);
	sys.AddFunction(x*y);

	Vec<complex_dbl> values(2);
	values << complex_dbl(2.0), complex_dbl(3.0);
	(void) sys.Eval(values);          // primes the SLP cache: [y+1, x*y]

	sys += sys;                       // mutates tree to [2(y+1), 2xy] but cache stale

	auto after = sys.Eval(values);    // expected [8, 12], currently [4, 6]
	BOOST_CHECK_EQUAL(after(0), complex_dbl(8.0));
	BOOST_CHECK_EQUAL(after(1), complex_dbl(12.0));
}


/**
\class bertini::System
\test \b operator_mult_equals_invalidates_slp_cache Sibling of the
operator+= cache-invalidation test. operator*= multiplies each function
by a Node and must also invalidate is_differentiated_.
*/
BOOST_AUTO_TEST_CASE(operator_mult_equals_invalidates_slp_cache)
{
	using bertini::node::Complex;

	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys.AddVariableGroup(vars);
	sys.AddFunction(x+y);

	Vec<complex_dbl> values(2);
	values << complex_dbl(1.0), complex_dbl(2.0);
	(void) sys.Eval(values);          // primes SLP for f = x+y

	sys *= Complex::Make("3.0");        // f should now be 3*(x+y)

	auto after = sys.Eval(values);    // expected [9], currently [3]
	BOOST_CHECK_EQUAL(after(0), complex_dbl(9.0));
}


/**
\class bertini::System
\test \b reorder_functions_decreasing_invalidates_slp_cache The function
ordering inside the SLP corresponds to the order of functions_ at SLP
build time. ReorderFunctionsByDegreeDecreasing swaps entries in functions_
and must invalidate is_differentiated_ so the SLP is rebuilt to match.
*/
BOOST_AUTO_TEST_CASE(reorder_functions_decreasing_invalidates_slp_cache)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);
	sys.AddVariableGroup(vars);

	sys.AddFunction(x+y);             // degree 1, initially at position 0
	sys.AddFunction(x*y);             // degree 2, initially at position 1

	Vec<complex_dbl> values(2);
	values << complex_dbl(2.0), complex_dbl(3.0);
	(void) sys.Eval(values);          // primes SLP with the [degree1, degree2] order

	sys.ReorderFunctionsByDegreeDecreasing();   // now [degree2, degree1] = [x*y, x+y]

	auto after = sys.Eval(values);
	BOOST_CHECK_EQUAL(after(0), complex_dbl(6.0));      // x*y at position 0
	BOOST_CHECK_EQUAL(after(1), complex_dbl(5.0));      // x+y at position 1
}


/**
\class bertini::System
\test \b reorder_functions_increasing_invalidates_slp_cache Sibling of the
decreasing test. Reorders an already-cached system into ascending degree
order and verifies the SLP follows.
*/
BOOST_AUTO_TEST_CASE(reorder_functions_increasing_invalidates_slp_cache)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);
	sys.AddVariableGroup(vars);

	sys.AddFunction(x*y);             // degree 2, initially at position 0
	sys.AddFunction(x+y);             // degree 1, initially at position 1

	Vec<complex_dbl> values(2);
	values << complex_dbl(2.0), complex_dbl(3.0);
	(void) sys.Eval(values);          // primes SLP with the [degree2, degree1] order

	sys.ReorderFunctionsByDegreeIncreasing();   // now [degree1, degree2] = [x+y, x*y]

	auto after = sys.Eval(values);
	BOOST_CHECK_EQUAL(after(0), complex_dbl(5.0));      // x+y at position 0
	BOOST_CHECK_EQUAL(after(1), complex_dbl(6.0));      // x*y at position 1
}


/**
\class bertini::System
\test \b eval_wrong_size_input_throws Verifies that passing a variable
vector of the wrong size to System::Eval throws std::runtime_error rather
than reading past the end or producing garbage.
*/
BOOST_AUTO_TEST_CASE(eval_wrong_size_input_throws)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys.AddVariableGroup(vars);
	sys.AddFunction(x+y);

	Vec<complex_dbl> too_small(1);
	too_small << complex_dbl(1.0);
	BOOST_CHECK_THROW(sys.Eval(too_small), std::runtime_error);

	Vec<complex_dbl> too_big(3);
	too_big << complex_dbl(1.0), complex_dbl(2.0), complex_dbl(3.0);
	BOOST_CHECK_THROW(sys.Eval(too_big), std::runtime_error);
}


/**
\class bertini::System
\test \b add_systems_chain Verifies operator+= return-by-reference and that
chained += accumulates correctly. C++ groups `a += b += c` as `a += (b += c)`,
so b is mutated to b+c, then a becomes a+b+c.
*/
BOOST_AUTO_TEST_CASE(add_systems_chain)
{
	bertini::System sys1, sys2, sys3;
	Var x = Variable::Make("x"), y = Variable::Make("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	for (auto* s : {&sys1, &sys2, &sys3})
		s->AddVariableGroup(vars);

	sys1.AddFunction(x);       sys1.AddFunction(y);
	sys2.AddFunction(y);       sys2.AddFunction(x);
	sys3.AddFunction(x*y);     sys3.AddFunction(x+y);

	sys1 += sys2 += sys3;      // sys2 becomes sys2+sys3; sys1 becomes sys1+sys2+sys3

	Vec<complex_dbl> v(2);
	v << complex_dbl(2.0), complex_dbl(3.0);

	// sys1 originally: [x, y] = [2, 3]
	// sys2 originally: [y, x] = [3, 2]
	// sys3 originally: [xy, x+y] = [6, 5]
	// sys1 final:      [2+3+6, 3+2+5] = [11, 10]
	auto v1 = sys1.Eval(v);
	BOOST_CHECK_EQUAL(v1(0), complex_dbl(11.0));
	BOOST_CHECK_EQUAL(v1(1), complex_dbl(10.0));

	// sys2 final:      [3+6, 2+5] = [9, 7]
	auto v2 = sys2.Eval(v);
	BOOST_CHECK_EQUAL(v2(0), complex_dbl(9.0));
	BOOST_CHECK_EQUAL(v2(1), complex_dbl(7.0));
}


/**
\class bertini::System
\test \b add_incompatible_systems_throws Exercises the four guard clauses in
System::operator+= for mismatched function count, variable count, and
variable group count. Each branch should throw std::runtime_error.
*/
BOOST_AUTO_TEST_CASE(add_incompatible_systems_throws)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars2; vars2.push_back(x); vars2.push_back(y);
	VariableGroup vars3; vars3.push_back(x); vars3.push_back(y); vars3.push_back(z);

	// Mismatched function count: 2 fns vs 1 fn, same vars.
	{
		bertini::System a, b;
		a.AddVariableGroup(vars2);  a.AddFunction(x);  a.AddFunction(y);
		b.AddVariableGroup(vars2);  b.AddFunction(x+y);
		BOOST_CHECK_THROW(a += b, std::runtime_error);
	}

	// Mismatched variable count: 2 vars vs 3 vars, same fn count.
	{
		bertini::System a, b;
		a.AddVariableGroup(vars2);  a.AddFunction(x+y);
		b.AddVariableGroup(vars3);  b.AddFunction(x+y+z);
		BOOST_CHECK_THROW(a += b, std::runtime_error);
	}

	// Mismatched variable group count: one VG of 2 vs two VGs of 1 each.
	{
		VariableGroup vg_x; vg_x.push_back(x);
		VariableGroup vg_y; vg_y.push_back(y);
		bertini::System a, b;
		a.AddVariableGroup(vars2);              a.AddFunction(x+y);
		b.AddVariableGroup(vg_x);  b.AddVariableGroup(vg_y);  b.AddFunction(x+y);
		BOOST_CHECK_THROW(a += b, std::runtime_error);
	}
}


/**
\class bertini::System
\test \b system_differentiate_wrt_time_linear Test the arithmetic sum of two Systems.
*/
BOOST_AUTO_TEST_CASE(system_differentiate_wrt_time_linear)
{
	Var x = Variable::Make("x");
	Var t = Variable::Make("t");
	auto f1 = (1-t)*x + t*(1-x);
	auto f2 = x-t;

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction(f1);
	S.AddFunction(f2);

	Vec<complex_dbl> v(1);
	v << 1.0;
	complex_dbl time(0.5,0.2);
	auto dS_dt = S.TimeDerivative(v,time);

	BOOST_CHECK_CLOSE( dS_dt(0).real(), complex_dbl(-1).real(), threshold_clearance_d);
	BOOST_CHECK_CLOSE( dS_dt(1).imag(), complex_dbl(-1).imag(), threshold_clearance_d);


}





/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_aff_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_aff_group)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x, y};
	sys.AddVariableGroup(vars);

	sys.Homogenize();

	Vec<complex_dbl> v(3);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),2);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_dehomogenize_FIFO_two_aff_groups Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_two_aff_groups)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	Var z = Variable::Make("z"), w = Variable::Make("w");
	VariableGroup vars{x, y};
	VariableGroup vars2{z, w};
	sys.AddVariableGroup(vars);
	sys.AddVariableGroup(vars2);

	sys.Homogenize();

	Vec<complex_dbl> v(6);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5), 
		 complex_dbl(5,6), complex_dbl(6,7), complex_dbl(7,8);

	auto d = sys.DehomogenizePoint(v);


	BOOST_CHECK_EQUAL(d.size(),4);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(4)/v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(5)/v(3)) < threshold_clearance_d);
}



/**
\class bertini::System
\test \b system_dehomogenize_FIFO_two_aff_groups_one_hom_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_two_aff_groups_one_hom_group)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	Var z = Variable::Make("z"), w = Variable::Make("w");
	Var h1 = Variable::Make("h1"), h2 = Variable::Make("h2");
	VariableGroup vars{x, y};
	VariableGroup vars2{h1,h2};
	VariableGroup vars3{z, w};
	sys.AddVariableGroup(vars);
	sys.AddHomVariableGroup(vars2);
	sys.AddVariableGroup(vars3);

	sys.Homogenize();

	Vec<complex_dbl> v(8);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5), 
		 complex_dbl(10,11), complex_dbl(11,12),
		 complex_dbl(5,6), complex_dbl(6,7), complex_dbl(7,8);

	auto d = sys.DehomogenizePoint(v);


	BOOST_CHECK_EQUAL(d.size(),6);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(4)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(4) - v(6)/v(5)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(5) - v(7)/v(5)) < threshold_clearance_d);
}



/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_hom_group Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_hom_group)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x, y};
	sys.AddHomVariableGroup(vars);

	sys.Homogenize();

	Vec<complex_dbl> v(2);
	v << complex_dbl(2,3), complex_dbl(3,4);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),2);

	BOOST_CHECK(abs(d(0) - v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(1)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_homogenize_point_unhomogenized_is_identity HomogenizePoint on an unhomogenized, unpatched system returns the point unchanged.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_point_unhomogenized_is_identity)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x, y};
	sys.AddVariableGroup(vars);

	Vec<complex_dbl> p(2);
	p << complex_dbl(3,4), complex_dbl(4,5);

	auto h = sys.HomogenizePoint(p);

	BOOST_CHECK_EQUAL(h.size(),2);
	BOOST_CHECK(abs(h(0) - p(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(1) - p(1)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_homogenize_point_FIFO_one_aff_group HomogenizePoint inserts the homogenizing coordinate with value 1, and dehomogenization inverts it.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_point_FIFO_one_aff_group)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x, y};
	sys.AddVariableGroup(vars);

	sys.Homogenize();

	Vec<complex_dbl> p(2);
	p << complex_dbl(3,4), complex_dbl(4,5);

	auto h = sys.HomogenizePoint(p);

	BOOST_CHECK_EQUAL(h.size(),3);
	BOOST_CHECK(abs(h(0) - complex_dbl(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(1) - p(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(2) - p(1)) < threshold_clearance_d);

	auto p_again = sys.DehomogenizePoint(h);
	BOOST_CHECK_EQUAL(p_again.size(),2);
	BOOST_CHECK(abs(p_again(0) - p(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(p_again(1) - p(1)) < threshold_clearance_d);
}


/**
\class bertini::System
\test \b system_homogenize_point_FIFO_mixed_groups HomogenizePoint handles affine groups, hom groups, and dehomogenization inverts it, with multiple group types in play.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_point_FIFO_mixed_groups)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	Var z = Variable::Make("z"), w = Variable::Make("w");
	Var h1 = Variable::Make("h1"), h2 = Variable::Make("h2");
	VariableGroup vars{x, y};
	VariableGroup vars2{h1,h2};
	VariableGroup vars3{z, w};
	sys.AddVariableGroup(vars);
	sys.AddHomVariableGroup(vars2);
	sys.AddVariableGroup(vars3);

	sys.Homogenize();

	Vec<complex_dbl> p(6);
	p << complex_dbl(3,4), complex_dbl(4,5),
		 complex_dbl(10,11), complex_dbl(11,12),
		 complex_dbl(6,7), complex_dbl(7,8);

	auto h = sys.HomogenizePoint(p);

	BOOST_CHECK_EQUAL(h.size(),8);
	// [hom0, x, y, hom-group passthrough, hom1, z, w]
	BOOST_CHECK(abs(h(0) - complex_dbl(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(1) - p(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(2) - p(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(3) - p(2)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(4) - p(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(5) - complex_dbl(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(6) - p(4)) < threshold_clearance_d);
	BOOST_CHECK(abs(h(7) - p(5)) < threshold_clearance_d);

	auto p_again = sys.DehomogenizePoint(h);
	BOOST_CHECK_EQUAL(p_again.size(),6);
	for (int ii = 0; ii < 6; ++ii)
		BOOST_CHECK(abs(p_again(ii) - p(ii)) < threshold_clearance_d);
}




/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_hom_group_two_ungrouped_vars Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_hom_group_two_ungrouped_vars)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	Var z = Variable::Make("z"), w = Variable::Make("w");
	VariableGroup vars{x, y};

	sys.AddHomVariableGroup(vars);
	sys.AddUngroupedVariable(z);
	sys.AddUngroupedVariable(w);

	sys.Homogenize();

	Vec<complex_dbl> v(4);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5), complex_dbl(5,6);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),4);

	BOOST_CHECK(abs(d(0) - v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(1)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(2) - v(2)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(3)) < threshold_clearance_d);
}

/**
\class bertini::System
\test \b system_dehomogenize_FIFO_one_aff_group_two_ungrouped_vars_another_aff_grp_hom_grp Test the dehomogenization of a point using the first-in-first-out variable ordering which is standard in Bertini 1.
*/
BOOST_AUTO_TEST_CASE(system_dehomogenize_FIFO_one_aff_group_two_ungrouped_vars_another_aff_grp_hom_grp)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	Var z = Variable::Make("z"), w = Variable::Make("w");
	Var h1 = Variable::Make("h1"), h2 = Variable::Make("h2");
	Var u1 = Variable::Make("u1"), u2 = Variable::Make("u2");

	VariableGroup vars{x,y};
	VariableGroup vars2{z,w};

	VariableGroup vars3{h1,h2};

	sys.AddVariableGroup(vars);
	sys.AddUngroupedVariable(u1);
	sys.AddUngroupedVariable(u2);
	sys.AddVariableGroup(vars2);
	sys.AddHomVariableGroup(vars3);

	sys.Homogenize();

	Vec<complex_dbl> v(10);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5), 
		 complex_dbl(10,11), complex_dbl(11,12),
		 complex_dbl(5,6), complex_dbl(6,7), complex_dbl(7,8),
		 complex_dbl(12,13), complex_dbl(13,14);

	auto d = sys.DehomogenizePoint(v);

	BOOST_CHECK_EQUAL(d.size(),8);

	BOOST_CHECK(abs(d(0) - v(1)/v(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(1) - v(2)/v(0)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(2) - v(3)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(3) - v(4)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(4) - v(6)/v(5)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(5) - v(7)/v(5)) < threshold_clearance_d);

	BOOST_CHECK(abs(d(6) - v(8)) < threshold_clearance_d);
	BOOST_CHECK(abs(d(7) - v(9)) < threshold_clearance_d);
}










/**
\class bertini::System
\test \b system_estimate_coeff_bound_linear Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_linear)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Var x = Variable::Make("x");
	Var t = Variable::Make("t");

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction((1-t)*x + t*(1-x));
	S.AddFunction(x-t);

	real_mp coefficient_bound = S.CoefficientBound<mpfr>();
	BOOST_CHECK(coefficient_bound < real_mp("10"));
	BOOST_CHECK(coefficient_bound > real_mp("0.5"));
}


/**
\class bertini::System
\test \b system_estimate_coeff_bound_quartic Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_quartic)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + real_mp("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	real_mp coefficient_bound = sys.CoefficientBound<mpfr>();
	BOOST_CHECK(coefficient_bound < real_mp("5"));
	BOOST_CHECK(coefficient_bound > real_mp("2"));
}



/**
\class bertini::System
\test \b system_estimate_coeff_bound_quartic Test the estimation of the largest coefficient in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_coeff_bound_homogenized_quartic)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	// AutoPatch coefficients come from RandomMp; seeding makes this test's draws
	// (and hence the <10 threshold below) deterministic regardless of which other
	// tests ran first.  RandomMp now honors SetGlobalSeed (it shares ThreadEngine).
	bertini::SetGlobalSeed(1u);

	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + real_mp("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	sys.Homogenize();
	sys.AutoPatch();

	real_mp coefficient_bound = sys.CoefficientBound<mpfr>();
	BOOST_CHECK(coefficient_bound < real_mp("10"));
	BOOST_CHECK(coefficient_bound > real_mp("2"));
}


/**
\class bertini::System
\test \b system_homogenize_point_lands_on_patch On a patched system, HomogenizePoint produces a point ON the patch (rescaling it again is the identity), and dehomogenizing recovers the user point.
*/
BOOST_AUTO_TEST_CASE(system_homogenize_point_lands_on_patch)
{
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x, y};
	sys.AddVariableGroup(vars);

	sys.Homogenize();
	sys.AutoPatch();

	Vec<complex_dbl> p(2);
	p << complex_dbl(3,4), complex_dbl(4,5);

	auto h = sys.HomogenizePoint(p);
	BOOST_CHECK_EQUAL(h.size(),3);

	// already on the patch: rescaling is the identity
	auto h_rescaled = sys.RescalePointToFitPatch(h);
	for (int ii = 0; ii < 3; ++ii)
		BOOST_CHECK(abs(h_rescaled(ii) - h(ii)) < threshold_clearance_d);

	// projectively the same point: dehomogenizing recovers the user coordinates
	auto p_again = sys.DehomogenizePoint(h);
	BOOST_CHECK(abs(p_again(0) - p(0)) < threshold_clearance_d);
	BOOST_CHECK(abs(p_again(1) - p(1)) < threshold_clearance_d);

	// and the round trip from an on-patch internal point is exact:
	// Hom(Dehom(s)) == s for s on the patch
	Vec<complex_dbl> v(3);
	v << complex_dbl(2,3), complex_dbl(3,4), complex_dbl(4,5);
	auto s = sys.RescalePointToFitPatch(v);
	auto s_again = sys.HomogenizePoint(sys.DehomogenizePoint(s));
	for (int ii = 0; ii < 3; ++ii)
		BOOST_CHECK(abs(s_again(ii) - s(ii)) < threshold_clearance_d);
}

/**
\class bertini::System
\test \b system_estimate_degree_bound_linear Test the estimation of the degree in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_degree_bound_linear)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	Var x = Variable::Make("x");
	Var t = Variable::Make("t");

	bertini::System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction((1-t)*x + t*(1-x));
	S.AddFunction(x-t);

	real_mp degree_bound = S.DegreeBound();
	BOOST_CHECK(degree_bound == real_mp("1"));
}


/**
\class bertini::System
\test \b system_estimate_degree_bound_quartic Test the estimation of the degree in a system, including its derivatives.
*/
BOOST_AUTO_TEST_CASE(system_estimate_degree_bound_quartic)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + real_mp("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	real_mp degree_bound = sys.DegreeBound();
	BOOST_CHECK(degree_bound == real_mp("4"));
}




/**
\class bertini::System
\test \b system_multiply_by_node Ensure that multiplication of a system by a node doesn't affect other copies of a system
*/
BOOST_AUTO_TEST_CASE(system_multiply_by_node)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	bertini::System sys1, sys2;
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	sys1.AddVariableGroup(vars);  
	sys1.AddFunction(x);
	sys1.AddFunction(y);
	sys1.AddFunction(z);

	sys2.AddVariableGroup(vars);  
	sys2.AddFunction(y+x*y + real_mp("0.5"));
	sys2.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys2.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	Var t = Variable::Make("t");

	auto sys_copy1 = t*sys1;
	auto sys_copy2 = (1-t)*sys2;

	auto sys_copy3 = sys_copy1 + sys_copy2; // this line couples the two systems...  this coupling is total garbage!!!  this 'homotopy' should never be used.
}



/**
\class bertini::System
\test \b concatenate_two_systems Test that contactenation of two systems works correctly
*/
BOOST_AUTO_TEST_CASE(concatenate_two_systems)
{
	bertini::System sys1, sys2;
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	sys1.AddVariableGroup(vars);  
	sys1.AddFunction(x);
	sys1.AddFunction(y);
	sys1.AddFunction(z);

	sys2.AddVariableGroup(vars);  
	sys2.AddFunction(y+x*y + real_mp("0.5"));
	sys2.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys2.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);


	auto sys3 = Concatenate(sys1, sys2);

	BOOST_CHECK_EQUAL(sys3.NumNaturalFunctions(),6);
}

/**
\class bertini::System
\test \b concatenate_accepts_a_structured_block_system  Concatenate must append the functions of a
system whose rows live in a structured block (a linear-forms slice), not only PolynomialBlock rows.
Regression: it iterated sys2.Function(ii) -- which reads only the PolynomialBlock -- so a structured
operand was skipped (and null-deref/segfaulted when it had no polynomial block).
*/
BOOST_AUTO_TEST_CASE(concatenate_accepts_a_structured_block_system)
{
	Var x = Variable::Make("x"), y = Variable::Make("y");
	bertini::VariableGroup vars{x,y};

	bertini::System sys1;
	sys1.AddVariableGroup(vars);
	sys1.AddFunction(x*x + y*y - 1);          // a polynomial: f0 = x^2 + y^2 - 1

	// a system whose two functions live in a LinearFormsBlock (a slice): f = 2x+3y+1, x-y+4
	bertini::Mat<bertini::complex_mp> M(2,3);
	M << bertini::complex_mp(2), bertini::complex_mp(3),  bertini::complex_mp(1),
	     bertini::complex_mp(1), bertini::complex_mp(-1), bertini::complex_mp(4);
	auto sys2 = bertini::Slice::FromCoefficients(vars, M).AsSystem();

	auto sys3 = Concatenate(sys1, sys2);      // used to segfault on the structured operand
	BOOST_CHECK_EQUAL(sys3.NumNaturalFunctions(), 3);

	// at (x,y) = (1,1):  f0 = 1,  f1 = 6,  f2 = 4
	bertini::Vec<bertini::complex_dbl> p(2); p << bertini::complex_dbl(1), bertini::complex_dbl(1);
	auto v = sys3.Eval(p);
	BOOST_CHECK_EQUAL(v.size(), 3);
	BOOST_CHECK_CLOSE(v(0).real(), 1.0, 1e-11);
	BOOST_CHECK_CLOSE(v(1).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(v(2).real(), 4.0, 1e-11);
}

/**
\class bertini::System
\test \b concatenate_copies_patch_from_patched_operand When exactly one operand is patched,
Concatenate must give the result that operand's patch.  Regression: Concatenate used to call
sys1.CopyPatches(sys1) -- copying patches from itself, a no-op -- so the patch was silently
dropped and the result came out unpatched.
*/
BOOST_AUTO_TEST_CASE(concatenate_copies_patch_from_patched_operand)
{
	bertini::SetGlobalSeed(1u); // AutoPatch draws random patch coefficients

	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	bertini::System base;
	base.AddVariableGroup(vars);
	base.AddFunction(x*x + y*y - 1);
	base.Homogenize();

	// Two systems with identical variable ordering (clones share the node DAG), only one patched.
	auto sys1 = Clone(base);              // unpatched
	auto sys2 = Clone(base);
	sys2.AutoPatch();                     // patched

	BOOST_CHECK(!sys1.IsPatched());
	BOOST_CHECK(sys2.IsPatched());

	auto sys3 = Concatenate(sys1, sys2);
	BOOST_CHECK(sys3.IsPatched());
	BOOST_CHECK(sys3.GetPatch() == sys2.GetPatch());
}

/**
\class bertini::System
\test \b parsed_system_evaluates_correctly 
*/
BOOST_AUTO_TEST_CASE(parsed_system_evaluates_correctly)
{
	
	std::string str = "function f; variable_group x1, x2; y = x1*x2; f = y*y;";
	
	bertini::System sys;
	[[maybe_unused]] bool s = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	
	Vec<complex_dbl> values(2);
	
	values(0) = complex_dbl(2.0);
	values(1) = complex_dbl(3.0);
	
	Vec<complex_dbl> v(sys.NumNaturalFunctions());
	sys.EvalInPlace(v, values);
	
	BOOST_CHECK_EQUAL(v(0), 36.0);
	
	
	auto J = sys.Jacobian(values);
	
	double x1 = 2;
	double x2 = 3;
	
	BOOST_CHECK_EQUAL(J(0,0), 2*x1*x2*x2);
	BOOST_CHECK_EQUAL(J(0,1), x1*x1*2*x2);
}






BOOST_AUTO_TEST_CASE(variable_group_sizes_and_degrees_homvargrp)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");

	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};

	sys.AddHomVariableGroup(v1);
	sys.AddHomVariableGroup(v2);

	sys.AddFunction(x*y - 1);
	sys.AddFunction(pow(x,2) - 1);

	auto size_of_each_var_gp = sys.VariableGroupSizes(); 
	
	BOOST_CHECK_EQUAL(size_of_each_var_gp[0], 1);
	BOOST_CHECK_EQUAL(size_of_each_var_gp[1], 1);
	
	BOOST_CHECK(!sys.IsHomogeneous());
}

BOOST_AUTO_TEST_CASE(variable_group_sizes_and_degrees_affvargrp)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");

	System sys;

	VariableGroup v1{x};
	VariableGroup v2{y};

	sys.AddVariableGroup(v1);
	sys.AddVariableGroup(v2);

	sys.AddFunction(x*y - 1);
	sys.AddFunction(pow(x,2) - 1);

	auto size_of_each_var_gp = sys.VariableGroupSizes(); 
	
	BOOST_CHECK_EQUAL(size_of_each_var_gp[0], 1);
	BOOST_CHECK_EQUAL(size_of_each_var_gp[1], 1);
	

	sys.Homogenize();

	size_of_each_var_gp = sys.VariableGroupSizes(); 
	
	BOOST_CHECK_EQUAL(size_of_each_var_gp[0], 2);
	BOOST_CHECK_EQUAL(size_of_each_var_gp[1], 2);
}



BOOST_AUTO_TEST_CASE(clone_system_new_variables_evaluation)
{
	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	auto sys = bertini::system::Precon::GriewankOsborn();
	Vec<mpfr> x1(2), x2(2);
	x1(0) = bertini::multiprecision::RandomUnit(CLASS_TEST_MPFR_DEFAULT_DIGITS);
	x1(1) = bertini::multiprecision::RandomUnit(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	auto f = sys.Eval(x1);

	auto sys_clone = bertini::Clone(sys);

	x2(0) = mpfr{2};
	x2(1) = mpfr{3};

	auto f_clone = sys_clone.Eval(x2);

	BOOST_CHECK_EQUAL(mpfr(2.5), f_clone(0));
	BOOST_CHECK_EQUAL(mpfr(-1), f_clone(1));

	auto f_clone2 = sys_clone.Eval(x1);
	auto f2 = sys.Eval<mpfr>();

	BOOST_CHECK_EQUAL(f,f2);
	BOOST_CHECK_EQUAL(f_clone2,f2);
}



/**
\class bertini::node
\test \b gather_variables_alphabetical GatherVariables returns the distinct
variables of a set of functions, de-duplicated and sorted alphabetically by name.
*/
BOOST_AUTO_TEST_CASE(gather_variables_alphabetical)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");

	// note: declared out of alphabetical order, x used twice
	auto f1 = (pow(z,2) + y*x);
	auto f2 = (x - y);

	auto found = bertini::node::GatherVariables(std::vector<std::shared_ptr<bertini::node::Node>>{f1, f2});

	BOOST_CHECK_EQUAL(found.size(), 3);
	BOOST_CHECK_EQUAL(found[0]->name(), "x");
	BOOST_CHECK_EQUAL(found[1]->name(), "y");
	BOOST_CHECK_EQUAL(found[2]->name(), "z");
}


/**
\class bertini::System
\test \b system_construct_from_functions Constructing a System from a list of
functions auto-discovers the variables into a single affine variable group.
*/
BOOST_AUTO_TEST_CASE(system_construct_from_functions)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");

	auto f1 = (x*y*z);
	auto f2 = (x + y + z);

	bertini::System sys(std::vector<std::shared_ptr<bertini::node::Node>>{f1, f2});

	BOOST_CHECK_EQUAL(sys.NumTotalFunctions(), 2);
	BOOST_CHECK_EQUAL(sys.NumVariableGroups(), 1);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 3);

	auto const& ordering = sys.Variables();
	BOOST_CHECK_EQUAL(ordering[0]->name(), "x");
	BOOST_CHECK_EQUAL(ordering[1]->name(), "y");
	BOOST_CHECK_EQUAL(ordering[2]->name(), "z");
}


/**
\class bertini::System
\test \b system_set_variable_groups SetVariableGroups replaces the whole
variable-group structure with the supplied affine groups.
*/
BOOST_AUTO_TEST_CASE(system_set_variable_groups)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var z = Variable::Make("z");

	auto f1 = (x*y*z);

	bertini::System sys(std::vector<std::shared_ptr<bertini::node::Node>>{f1});
	BOOST_CHECK_EQUAL(sys.NumVariableGroups(), 1);

	bertini::VariableGroup g1{x};
	bertini::VariableGroup g2{y, z};
	sys.SetVariableGroups(std::vector<bertini::VariableGroup>{g1, g2});

	BOOST_CHECK_EQUAL(sys.NumVariableGroups(), 2);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 3);
}


BOOST_AUTO_TEST_SUITE_END()




