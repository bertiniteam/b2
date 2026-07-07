//This file is part of Bertini 2.
//
//b2/core/test/classes/function_tree_transform.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/core/test/classes/function_tree_transform.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/core/test/classes/function_tree_transform.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  silviana amethyst, university of wisconsin-eau claire

/**
\file Testing suite for transforms applied to the function tree.  In this suite, we content ourselves to evaluation in double precision, exercising multiple precision evaluation elsewhere.
*/

#include <iostream>
#include <sstream>

#include <cstdlib>
#include <cmath>

#include "bertini2/function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include "externs.hpp"
#include "eval_helper.hpp"

using bertini::test::EvalAt;

using Nd = std::shared_ptr<bertini::node::Node>;

using Variable = bertini::node::Variable;
using Node = bertini::node::Node;
using Integer = bertini::node::Integer;
using Rational = bertini::node::Rational;
using SumOperator = bertini::node::SumOperator;
using MultOperator = bertini::node::MultOperator;

using complex_dbl = bertini::complex_dbl;

auto MakeZero(){return Nd(Integer::Make(0));}
auto MakeOne(){return Nd(Integer::Make(1));}

// the printed form of a simplified expression
inline std::string SimplifiedForm(Nd const& n){ std::ostringstream o; bertini::Simplify(n)->print(o); return o.str(); }

BOOST_AUTO_TEST_SUITE(function_tree)

BOOST_AUTO_TEST_SUITE(transform)


//       _ _           _             _                               
//      | (_)         (_)           | |                              
//   ___| |_ _ __ ___  _ _ __   __ _| |_ ___   _______ _ __ ___  ___ 
//  / _ \ | | '_ ` _ \| | '_ \ / _` | __/ _ \ |_  / _ \ '__/ _ \/ __|
// |  __/ | | | | | | | | | | | (_| | ||  __/  / /  __/ | | (_) \__ \.
//  \___|_|_|_| |_| |_|_|_| |_|\__,_|\__\___| /___\___|_|  \___/|___/
                                                                  
                                                                  


BOOST_AUTO_TEST_SUITE(simplified)

// Node::Simplified() is the non-mutating successor to the old in-place
// EliminateZeros / EliminateOnes / ReduceDepth machinery: it returns a NEW simplified
// tree and never touches its input.  Only *literal* zeros/ones are folded (per ADR-0011);
// an expression that merely happens to evaluate to zero is deliberately left alone.

// ---- literal zeros vanish / collapse ----
BOOST_AUTO_TEST_CASE(leaf_zero_is_unchanged)
{
	auto zero = MakeZero();
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(zero->Simplified()), 0.);
}

BOOST_AUTO_TEST_CASE(zeros_drop_from_sum_signs_preserved)
{
	auto zero = MakeZero();
	auto n = 2 + zero - 1 + zero - 2;
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(n)), -1.);
}

BOOST_AUTO_TEST_CASE(sums_of_zeros_are_zero)
{
	auto zero = MakeZero();
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(zero+zero)), 0.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(zero+0)), 0.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(0+zero)), 0.);
}

// ---- like factors combine into powers, like terms into coefficients ----
// Identical subexpressions are one interned node, so "structurally equal" is pointer identity:
// grouping factors/terms is a hash on the operand pointer.

BOOST_AUTO_TEST_CASE(like_factors_combine_into_powers)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_EQUAL(SimplifiedForm(x*x), "x^2");
	BOOST_CHECK_EQUAL(SimplifiedForm(x*x*x), "x^3");
	BOOST_CHECK_EQUAL(SimplifiedForm(pow(x,2)*pow(x,3)), "x^5");   // x^a * x^b -> x^(a+b)
	BOOST_CHECK_EQUAL(SimplifiedForm((x+y)*(x+y)), "(x+y)^2");     // any repeated base, not just vars
}

BOOST_AUTO_TEST_CASE(divided_like_factors_lower_the_exponent)
{
	auto x = Variable::Make("x");
	BOOST_CHECK_EQUAL(SimplifiedForm(x*x/x), "x");                 // x^2 / x -> x
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(x/x)), 1.);    // x / x -> 1
}

BOOST_AUTO_TEST_CASE(like_terms_combine_into_coefficients)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_EQUAL(SimplifiedForm(x+x), "2*x");
	BOOST_CHECK_EQUAL(SimplifiedForm(x+x+x), "3*x");
	BOOST_CHECK_EQUAL(SimplifiedForm(Integer::Make(3)*x + Integer::Make(2)*x), "5*x");
	BOOST_CHECK_EQUAL(SimplifiedForm(x*y + x*y), "2*x*y");
}

BOOST_AUTO_TEST_CASE(opposite_like_terms_cancel)
{
	auto x = Variable::Make("x");
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(x - x)), 0.);
	BOOST_CHECK_EQUAL(SimplifiedForm(Integer::Make(3)*x - x), "2*x");
}

BOOST_AUTO_TEST_CASE(combining_preserves_value)
{
	auto x = Variable::Make("x");
	Nd f = x*x*x + x*x*x;                 // 2*x^3 = 16 at x=2
	auto s = bertini::Simplify(f);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(s, {{"x", complex_dbl(2.0, 0.0)}}), complex_dbl(16.0, 0.0));
}

BOOST_AUTO_TEST_CASE(zero_terms_drop_keeping_value)
{
	auto x = Variable::Make("x");
	auto zero = MakeZero();
	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(2.0)} };
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(pow(x,2) + zero*2), pt), 4.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(pow(x,2) + x*zero*2*(x*x*x)), pt), 4.);
}

BOOST_AUTO_TEST_CASE(zero_term_drops_in_deeper_sum)
{
	auto x = Variable::Make("x");
	auto zero = MakeZero();
	auto n = (x + sqrt(x) + zero) + x*2*(x*x*x);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(n), {{"x", complex_dbl(2.0)}}),
		35.414213562373095048801688724209698078569671875377);
}

BOOST_AUTO_TEST_CASE(literal_zero_collapses_product)
{
	auto zero = MakeZero();
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(1*zero)), 0.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(zero*zero)), 0.);
}

// ---- literal ones drop from products ----
BOOST_AUTO_TEST_CASE(leaf_one_is_unchanged)
{
	auto one = MakeOne();
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(one->Simplified()), 1.);
}

BOOST_AUTO_TEST_CASE(ones_drop_from_product)
{
	auto x = Variable::Make("x");
	auto one = MakeOne();
	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(2.)} };
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(one*one)), 1.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(x*one), pt), 2.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(one*x), pt), 2.);
}

BOOST_AUTO_TEST_CASE(ones_fold_through_division_chains)
{
	auto x = Variable::Make("x");
	auto one = Integer::Make(1);
	auto two = Integer::Make(2);
	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(2.)} };
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify((two*one/two)*x), pt), 2.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify(one*one/one*one/one*one*x), pt), 2.);
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(bertini::Simplify((one*one*2) * (one*x*one)), pt), 4.);
}

// ---- nested singletons unwrap to the variable itself ----
BOOST_AUTO_TEST_CASE(nested_singleton_mults_unwrap_to_the_variable)
{
	auto x = Variable::Make("x");
	Nd n = MultOperator::Make(x);
	for (int i = 0; i < 4; ++i) n = MultOperator::Make(n);
	auto s = bertini::Simplify(n);
	BOOST_CHECK(std::dynamic_pointer_cast<bertini::node::Variable>(s));  // unwrapped
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(s, {{"x", complex_dbl(5.)}}), 5.);
}

// ---- the load-bearing contract: Simplified() does NOT mutate its input ----
BOOST_AUTO_TEST_CASE(simplify_does_not_mutate_input)
{
	auto x = Variable::Make("x");
	auto zero = MakeZero();
	auto n = std::dynamic_pointer_cast<bertini::node::NaryOperator>((x + zero) + x);
	BOOST_REQUIRE(n);
	auto before = n->NumOperands();
	auto s = bertini::Simplify(n);
	BOOST_CHECK_EQUAL(n->NumOperands(), before);   // input untouched
	BOOST_CHECK(s != n);                            // a fresh tree was returned
	BOOST_CHECK_EQUAL(EvalAt<complex_dbl>(s, {{"x", complex_dbl(3.)}}), 6.);   // value preserved (x + x)
}

BOOST_AUTO_TEST_SUITE_END() // simplified














//        .__               .__  .__  _____       
//   _____|__| _____ ______ |  | |__|/ ____\__.__.
//  /  ___/  |/     \\____ \|  | |  \   __<   |  |
//  \___ \|  |  Y Y  \  |_> >  |_|  ||  |  \___  |
// /____  >__|__|_|  /   __/|____/__||__|  / ____|
//      \/         \/|__|                  \/     
//
//
// http://patorjk.com/software/taag/#p=display&f=Graffiti&t=Type%20Something%20







BOOST_AUTO_TEST_SUITE(simplify)

BOOST_AUTO_TEST_CASE(flattens_completely)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd m = SumOperator::Make(x, true);
	Nd n = SumOperator::Make(y, true);

	auto p = m+n+0;
	auto q = 0+m-n*1;
	auto s = pow(x,2);

	auto r = p+q+0*s + 0*1 + sqrt(0*x);

complex_dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
complex_dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	auto rs = bertini::Simplify(r);  // functional: r is untouched, rs is simplified

	auto result = EvalAt<complex_dbl>(rs, {{"x", a}, {"y", b}});
	BOOST_CHECK_EQUAL(result, a+a);
}

BOOST_AUTO_TEST_CASE(complicated)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	auto n = (((((2*x*1)*y)+(0*(pow(x,2))))/2)-(0*((pow(x,2))*y)/(static_cast<int>(pow(2,2)))));

	auto ns = bertini::Simplify(n);

	complex_dbl a(1.7, -0.3);
	complex_dbl b(-0.9, 2.1);
	std::map<std::string,complex_dbl> pt{ {"x", a}, {"y", b} };

	BOOST_CHECK_SMALL(abs(EvalAt<complex_dbl>(ns, pt) - a*b), threshold_clearance_d);
}


BOOST_AUTO_TEST_CASE(complicated2)
{
	auto x = Variable::Make("x");
	auto t = Variable::Make("t");

	auto f = (((pow((x-1),2))*(1-t))+((pow(x,2)+1)*t));

	auto dfdx = f->Differentiate(x);

	auto dfdxs = bertini::Simplify(dfdx);

	complex_dbl xval(1.3, -0.7);
	complex_dbl tval(0.4, 0.2);
	std::map<std::string,complex_dbl> pt{ {"x", xval}, {"t", tval} };

	BOOST_CHECK_SMALL(abs(EvalAt<complex_dbl>(dfdxs, pt) - 2.*(xval+tval-1.)), 1e-15);
}


BOOST_AUTO_TEST_CASE(complicated3)
{
	auto zero = MakeZero();
	auto one = MakeOne();

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto HOM_VAR_0 = Variable::Make("HOM_VAR_0");
	auto t = Variable::Make("t");

	auto f = ((zero*pow((y-(HOM_VAR_0*one)),2))+((2*(y-(HOM_VAR_0*one))*(one-((zero*one)+(zero*HOM_VAR_0))))*(one-t)));

	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(1.1,-0.2)}, {"y", complex_dbl(0.7,1.3)},
		{"HOM_VAR_0", complex_dbl(-0.5,0.9)}, {"t", complex_dbl(0.3,0.6)} };

	auto init_val = EvalAt<complex_dbl>(f, pt);

	auto fs = bertini::Simplify(f);

	// simplify reorders the arithmetic (like-factor/like-term combining), so compare with a
	// tolerance rather than for bit-exact equality.
	BOOST_CHECK_SMALL(std::abs(init_val - EvalAt<complex_dbl>(fs, pt)), 1e-12);

}

BOOST_AUTO_TEST_CASE(complicated4)
{
	auto zero = MakeZero();
	auto one = MakeOne();

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto HOM_VAR_0 = Variable::Make("HOM_VAR_0");
	auto t = Variable::Make("t");

	complex_dbl a(1.4, -0.6);
	complex_dbl h(-0.8, 0.5);
	complex_dbl T(0.25, 0.7);
	std::map<std::string,complex_dbl> pt{ {"x", a}, {"y", complex_dbl(0.3,0.1)}, {"HOM_VAR_0", h}, {"t", T} };

	auto f = ((zero*(pow((x-(HOM_VAR_0*one)),3)))+((3*(pow((x-(HOM_VAR_0*one)),2))*(-((one*one)+(zero*HOM_VAR_0))))*(one-t)));

	auto f_val_init = EvalAt<complex_dbl>(f, pt);
	[[maybe_unused]] auto actual_val = 3.*pow((h - a),2)*(T - 1.);

	auto fs = bertini::Simplify(f);

	auto f_val_after = EvalAt<complex_dbl>(fs, pt);
	// canonical operand ordering reorders sums/products, so simplify preserves the
	// value only up to reorder rounding -- compare with a tolerance, not exact equality.
	BOOST_CHECK_SMALL(std::abs(f_val_init - f_val_after), 1e-12);
}

//((0*((x-(HOM_VAR_0*1))^3))+((3*((x-(HOM_VAR_0*1))^2)*(-((1*1)+(0*HOM_VAR_0))))*(1-t)))

BOOST_AUTO_TEST_CASE(yet_more_complicated)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto t = Variable::Make("t");
	auto HOM_VAR_0 = Variable::Make("HOM_VAR_0");

	auto f = (((0*(pow(y,2)-(pow(HOM_VAR_0,2)*Rational::Make("9215126146405988386300422552813438491596469014004/18831809439874092151531390861220941995712612447011","-34979570316540871529966189550041755413443953059471/3298415588464117388268211317293113094514010909093"))))+(((2*y*1)-(((2*HOM_VAR_0*0)*Rational::Make("9215126146405988386300422552813438491596469014004/18831809439874092151531390861220941995712612447011","-34979570316540871529966189550041755413443953059471/3298415588464117388268211317293113094514010909093"))+(0*pow(HOM_VAR_0,2))))*t))+((0*pow((y-(HOM_VAR_0*1)),2))+((2*(y-(HOM_VAR_0*1))*(1-((0*1)+(0*HOM_VAR_0))))*(1-t))));

	// jac_fn(1,2) = (((0*((y^2)-((HOM_VAR_0^2)*(9215126146405988386300422552813438491596469014004/18831809439874092151531390861220941995712612447011,-34979570316540871529966189550041755413443953059471/3298415588464117388268211317293113094514010909093))))+(((2*y*1)-(((2*HOM_VAR_0*0)*(9215126146405988386300422552813438491596469014004/18831809439874092151531390861220941995712612447011,-34979570316540871529966189550041755413443953059471/3298415588464117388268211317293113094514010909093))+(0*(HOM_VAR_0^2))))*t))+((0*((y-(HOM_VAR_0*1))^2))+((2*(y-(HOM_VAR_0*1))*(1-((0*1)+(0*HOM_VAR_0))))*(1-t))))

	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(0.6,-0.4)}, {"y", complex_dbl(1.2,0.8)},
		{"t", complex_dbl(0.35,0.15)}, {"HOM_VAR_0", complex_dbl(-0.7,0.3)} };

	auto init_val = EvalAt<complex_dbl>(f, pt);

	auto fs = bertini::Simplify(f);
	// simplify combines like factors/terms, which reorders the arithmetic, so it preserves the
	// value only up to floating-point reorder rounding -- compare with a tolerance.
	BOOST_CHECK_SMALL(std::abs(init_val - EvalAt<complex_dbl>(fs, pt)), 1e-12);
}

BOOST_AUTO_TEST_SUITE_END() // simplify



BOOST_AUTO_TEST_SUITE_END() // transform
BOOST_AUTO_TEST_SUITE_END() // function_tree






