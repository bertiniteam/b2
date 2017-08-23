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
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  Dani Brake
//  University of Wisconsin Eau Claire
//  ACMS
//  Fall 2017

/**
\file Testing suite for transforms applied to the function tree.  In this suite, we content ourselves to evaluation in double precision, exercising multiple precision evaluation elsewhere.
*/

#include <iostream>

#include <cstdlib>
#include <cmath>

#include "bertini2/function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include "externs.hpp"

using Nd = std::shared_ptr<bertini::node::Node>;

using bertini::MakeVariable;
using bertini::MakeInteger;

using dbl = bertini::dbl;

auto MakeZero(){return Nd(MakeInteger(0));}
auto MakeOne(){return Nd(MakeInteger(1));}

BOOST_AUTO_TEST_SUITE(function_tree)

BOOST_AUTO_TEST_SUITE(transform)

BOOST_AUTO_TEST_SUITE(eliminate_zeros)

BOOST_AUTO_TEST_SUITE(sum)

BOOST_AUTO_TEST_CASE(cant_eliminate_self)
{
	auto zero = MakeZero();
	auto num_eliminated = zero->EliminateZeros();
	BOOST_CHECK_EQUAL(num_eliminated, 0);
	BOOST_CHECK_EQUAL(zero->Eval<dbl>(), 0.);
}

BOOST_AUTO_TEST_CASE(cant_eliminate_self2)
{
	auto zero = MakeZero();
	auto n = zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK_EQUAL(num_eliminated, 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}



BOOST_AUTO_TEST_CASE(level_one_signs_preserved)
{
	auto zero = MakeZero();
	auto n = 2 + zero - 1 + zero - 2;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), -1.);
}


BOOST_AUTO_TEST_CASE(level_one)
{
	auto zero = MakeZero();
	auto n = zero+zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}


BOOST_AUTO_TEST_CASE(level_one2)
{
	auto zero = MakeZero();
	auto n = zero+0;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}


BOOST_AUTO_TEST_CASE(level_one3)
{
	auto zero = MakeZero();
	auto n = 0+zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}






BOOST_AUTO_TEST_CASE(level_one_more_complicated)
{
	auto x = MakeVariable("x");
	auto zero = MakeZero();


	auto n = pow(x,2) + zero*2;

	x->set_current_value(dbl(2.0));

	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 4.);
}



BOOST_AUTO_TEST_CASE(level_one_more_complicated2)
{
	auto x = MakeVariable("x");
	auto zero = MakeZero();


	auto n = pow(x,2) + x*zero*2*(x*x*x);

	x->set_current_value(dbl(2.0));

	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 4.);
}

BOOST_AUTO_TEST_CASE(level_two)
{
	auto x = MakeVariable("x");
	auto zero = MakeZero();


	auto n = (x + sqrt(x) + zero) + x*2*(x*x*x);

	x->set_current_value(dbl(2.0));

	auto num_eliminated = n->EliminateZeros();

	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 35.414213562373095048801688724209698078569671875377);
}

BOOST_AUTO_TEST_CASE(level_two_variable_set_to_zero_eliminated)
{
	auto x = MakeVariable("x");
	auto zero = MakeZero();


	auto n = (x + sqrt(x) + zero) + x*2*(x*x*x);

	auto as_op = std::dynamic_pointer_cast<bertini::node::NaryOperator>(n);


	x->set_current_value(dbl(0));

	auto num_eliminated = n->EliminateZeros();

	x->set_current_value(dbl(2.0));

	BOOST_CHECK_EQUAL(as_op->children_size(), 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.0);
}


BOOST_AUTO_TEST_SUITE_END() // eliminate zeros :: sum



/////////
//
//    ppppppppp
//     pp   pp
//     pp   pp
//     pp   pp
//     pp   pp
//
//////////////

BOOST_AUTO_TEST_SUITE(product)


BOOST_AUTO_TEST_CASE(level_one)
{
	auto zero = MakeZero();
	auto n = 1*zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}


BOOST_AUTO_TEST_CASE(level_one2)
{
	auto zero = MakeZero();
	auto n = 1*zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}

BOOST_AUTO_TEST_CASE(level_one3)
{
	auto zero = MakeZero();
	auto n = zero*zero;
	auto num_eliminated = n->EliminateZeros();
	BOOST_CHECK(num_eliminated >= 1);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 0.);
}

BOOST_AUTO_TEST_SUITE_END() // eliminate zeros :: product


BOOST_AUTO_TEST_SUITE_END() // eliminate zeros











//////////////////////////////
//
//   ONES
//
//      1
//     11
//    1 1
//      1
//      1
//    11111
//
/////////////////////////

BOOST_AUTO_TEST_SUITE(eliminate_ones)
BOOST_AUTO_TEST_SUITE(product)

BOOST_AUTO_TEST_CASE(cant_eliminate_self)
{
	auto one = MakeOne();
	auto num_eliminated = one->EliminateOnes();
	BOOST_CHECK_EQUAL(num_eliminated, 0);
	BOOST_CHECK_EQUAL(one->Eval<dbl>(), 1.);
}

BOOST_AUTO_TEST_CASE(cant_eliminate_self2)
{
	auto one = MakeOne();
	auto n = one;
	auto num_eliminated = n->EliminateOnes();
	BOOST_CHECK_EQUAL(num_eliminated, 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 1.);
}




BOOST_AUTO_TEST_CASE(level_one)
{
	auto one = MakeOne();
	auto n = one*one;
	auto num_eliminated = n->EliminateOnes();
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 1.);
}


BOOST_AUTO_TEST_CASE(level_one2)
{
	auto x = MakeVariable("x");
	auto one = MakeOne();
	auto n = x*one;
	auto num_eliminated = n->EliminateOnes();
	x->set_current_value(2.);
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 2.);
}


BOOST_AUTO_TEST_CASE(level_one3)
{
	auto x = MakeVariable("x");
	auto one = MakeOne();
	auto n = one*x;
	auto num_eliminated = n->EliminateOnes();
	x->set_current_value(2.);
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 2.);
}


BOOST_AUTO_TEST_CASE(level_one4)
{
	auto x = MakeVariable("x");
	auto one = MakeInteger(1);
	auto two = MakeInteger(2);
	auto n = (two*one/two)*x;
	auto num_eliminated = n->EliminateOnes();
	x->set_current_value(2.);
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 2.);
}


BOOST_AUTO_TEST_CASE(level_one5)
{
	auto x = MakeVariable("x");
	auto one = MakeInteger(1);
	auto n = one*one/one*one/one*one*x;
	auto num_eliminated = n->EliminateOnes();
	x->set_current_value(2.);
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 2.);
}


BOOST_AUTO_TEST_CASE(level_two)
{
	auto x = MakeVariable("x");
	auto one = MakeInteger(1);
	auto n = (one*one*2) * (one*x*one);
	x->set_current_value(2.);

	auto num_eliminated = n->EliminateOnes();
	
	BOOST_CHECK(num_eliminated > 0);
	BOOST_CHECK_EQUAL(n->Eval<dbl>(), 4.);
}


BOOST_AUTO_TEST_SUITE_END() // predict
BOOST_AUTO_TEST_SUITE_END() // eliminate ones




BOOST_AUTO_TEST_SUITE(reduce_depth)

BOOST_AUTO_TEST_SUITE(sum)

BOOST_AUTO_TEST_CASE(eliminates_sum_of_single)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::SumOperator>(x, true);
	Nd n = std::make_shared<bertini::node::SumOperator>(y, true);

	auto p = m+n;

	unsigned num_eliminated = p->ReduceDepth();

	BOOST_CHECK(num_eliminated > 0);

	dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
	dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	x->set_current_value(a); y->set_current_value(b);

	auto result = p->Eval<dbl>();
	BOOST_CHECK_EQUAL(result, a+b);
}


BOOST_AUTO_TEST_CASE(double_sum_signs_distribute)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::SumOperator>(x, true);
	Nd n = std::make_shared<bertini::node::SumOperator>(y, true);

	auto p = m+n;
	auto q = m-n;

	auto r = p+q;

	{
		unsigned num_eliminated = r->ReduceDepth();
		BOOST_CHECK(num_eliminated > 0);
		auto R = std::dynamic_pointer_cast<bertini::node::SumOperator>(r);
		BOOST_CHECK_EQUAL(R->children_size(), 4);

		dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
		dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

		x->set_current_value(a); y->set_current_value(b);
		auto result = r->Eval<dbl>();
		BOOST_CHECK_EQUAL(result, (a+b) + (a-b));
	}

	{
		unsigned num_eliminated = r->ReduceDepth();
		BOOST_CHECK(num_eliminated > 0);
		auto R = std::dynamic_pointer_cast<bertini::node::SumOperator>(r);
		BOOST_CHECK_EQUAL(R->children_size(), 4);

		dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
		dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

		x->set_current_value(a); y->set_current_value(b);
		r->Reset();
		auto result = r->Eval<dbl>();
		BOOST_CHECK_EQUAL(result, (a+b) + (a-b));
	}

}



BOOST_AUTO_TEST_CASE(eliminates_mult_of_single)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::MultOperator>(x);
	Nd n = std::make_shared<bertini::node::MultOperator>(y);

	auto p = m+n;

	unsigned num_eliminated = p->ReduceDepth();

	BOOST_CHECK(num_eliminated > 0);

	dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
	dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	x->set_current_value(a); y->set_current_value(b);

	auto result = p->Eval<dbl>();
	BOOST_CHECK_EQUAL(result, a+b);
}
BOOST_AUTO_TEST_SUITE_END() // sum


BOOST_AUTO_TEST_SUITE(prod)

BOOST_AUTO_TEST_CASE(eliminates_sum_of_single)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::SumOperator>(x, true);
	Nd n = std::make_shared<bertini::node::SumOperator>(y, true);

	auto p = m*n;

	unsigned num_eliminated = p->ReduceDepth();

	BOOST_CHECK(num_eliminated > 0);

	dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
	dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	x->set_current_value(a); y->set_current_value(b);

	auto result = p->Eval<dbl>();
	BOOST_CHECK_EQUAL(result, a*b);
}


BOOST_AUTO_TEST_CASE(eliminates_mult_of_single)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::MultOperator>(x);
	Nd n = std::make_shared<bertini::node::MultOperator>(y);

	auto p = m*n;

	unsigned num_eliminated = p->ReduceDepth();

	BOOST_CHECK(num_eliminated > 0);

	dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
	dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	x->set_current_value(a); y->set_current_value(b);

	auto result = p->Eval<dbl>();
	BOOST_CHECK_EQUAL(result, a*b);
}



BOOST_AUTO_TEST_CASE(distributive)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");


	auto m = x*y;
	auto n = x/y;

	auto p = m/n; // should reduce to y*y

	unsigned num_eliminated = p->ReduceDepth();

	BOOST_CHECK(num_eliminated > 0);

	dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
	dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);

	x->set_current_value(a); y->set_current_value(b);

	auto result = p->Eval<dbl>();
	BOOST_CHECK_SMALL(abs(result - b*b), 1e-13);
}



BOOST_AUTO_TEST_SUITE_END() // prod

BOOST_AUTO_TEST_SUITE_END() // reduce_depth


BOOST_AUTO_TEST_SUITE(simplify)

BOOST_AUTO_TEST_CASE(flattens_completely)
{
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");

	Nd m = std::make_shared<bertini::node::SumOperator>(x, true);
	Nd n = std::make_shared<bertini::node::SumOperator>(y, true);

	auto p = m+n+0;
	auto q = 0+m-n*1;
	auto s = pow(x,2);

	auto r = p+q+0*s + 0*1 + sqrt(0*x);

dbl a(4.1203847861962345182734, -5.1234768951256847623781614314);
dbl b(-8.98798649152356714919234, 0.49879892634876018735619234);
x->set_current_value(a); y->set_current_value(b);

	auto num_rounds = bertini::Simplify(r);

	BOOST_CHECK(num_rounds >= 2);


	r->Reset();
	auto result = r->Eval<dbl>();
	BOOST_CHECK_EQUAL(result, a+a);
}


BOOST_AUTO_TEST_SUITE_END() // simplify



BOOST_AUTO_TEST_SUITE_END() // transform
BOOST_AUTO_TEST_SUITE_END() // function_tree






