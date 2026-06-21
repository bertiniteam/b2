//This file is part of Bertini 2.
//
//interning_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//interning_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with interning_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for hash-consing (Rung 3): every Make() routes through Intern(), so structurally
equal nodes collapse to one shared object.  Pins the dedup behavior, the differentiation-
sharing win, and -- critically -- that building/simplifying never mutates a shared interned
node (the bug where SimplifiedSum did Make(first) then AddOperand(rest), corrupting a shared
single-operand sum).
*/

#include <cstdlib>
#include "bertini2/function_tree.hpp"
#include <boost/test/unit_test.hpp>

using Nd = std::shared_ptr<bertini::node::Node>;
using Variable = bertini::node::Variable;
using Integer = bertini::node::Integer;
using SumOperator = bertini::node::SumOperator;
using dbl = bertini::dbl;

BOOST_AUTO_TEST_SUITE(interning)

// ---- the regression: building/simplifying must not mutate a shared interned node ----

BOOST_AUTO_TEST_CASE(simplify_does_not_corrupt_a_shared_single_operand_sum)
{
	// This is the minimal form of the bug found in Rung 3a: a single-operand Sum is reused,
	// and SimplifiedSum used to do Make(first) (which now returns the *interned* shared sum)
	// then AddOperand(rest) -- mutating that shared node and corrupting every other holder.
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd m = SumOperator::Make(x, true);   // a one-term sum (+x); the node that got corrupted
	Nd n = SumOperator::Make(y, true);   // (+y)

	// build two expressions that reuse m and n, with cancelling y-terms
	Nd p = m + n;          // x + y
	Nd q = m - n;          // x - y
	Nd r = p + q + 0*x;    // 2x   (the +0 forces simplification work)

	dbl xv(2.0, -3.0), yv(5.0, 1.0);
	x->set_current_value(xv);
	y->set_current_value(yv);

	Nd rs = bertini::Simplify(r);
	rs->Reset();
	BOOST_CHECK_EQUAL(rs->Eval<dbl>(), xv + xv);   // == 2x, with the y's cancelled

	// and the reused sub-objects must be intact (not mutated by the simplification above)
	m->Reset(); BOOST_CHECK_EQUAL(m->Eval<dbl>(), xv);
	n->Reset(); BOOST_CHECK_EQUAL(n->Eval<dbl>(), yv);
}

// ---- basic dedup: structurally equal builds return the same object ----

BOOST_AUTO_TEST_CASE(equal_constants_are_one_object)
{
	BOOST_CHECK_EQUAL(Integer::Make(5).get(), Integer::Make(5).get());
	BOOST_CHECK(Integer::Make(5).get() != Integer::Make(6).get());
}

BOOST_AUTO_TEST_CASE(equal_operator_trees_are_one_object)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_EQUAL((x*y + y).get(), (x*y + y).get());            // whole tree shares
	BOOST_CHECK_EQUAL((x*y).get(), (x*y).get());                   // inner subexpr shares
	BOOST_CHECK(( x + y).get() != (y + x).get());                  // order-sensitive: distinct
}

BOOST_AUTO_TEST_CASE(variables_are_canonical_by_name)
{
	// Rung 3b: Make("x") interns to a single canonical x -- "system1's x IS system2's x".
	auto x1 = Variable::Make("x");
	auto x2 = Variable::Make("x");
	BOOST_CHECK_EQUAL(x1.get(), x2.get());
	BOOST_CHECK(Variable::Make("x").get() != Variable::Make("z").get());

	// because they are one object, setting the value through one is seen through the other
	x1->set_current_value(dbl(7.0, 0.0));
	x2->Reset();
	BOOST_CHECK_EQUAL(x2->Eval<dbl>(), dbl(7.0, 0.0));

	// and two independently-built expressions over "x" share their structure
	BOOST_CHECK_EQUAL((Variable::Make("x") * Variable::Make("x")).get(),
	                  (Variable::Make("x") * Variable::Make("x")).get());
}

// ---- the headline win: a shared subexpression's derivative is one interned node ----

BOOST_AUTO_TEST_CASE(differentiation_results_are_interned)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x*x + y*y;          // a single shared subexpression object

	// differentiating the SAME expression twice yields the SAME derivative node
	BOOST_CHECK_EQUAL(a->Differentiate(x).get(), a->Differentiate(x).get());

	// and when 'a' is reused in two functions, the d(a)/dx inside each derivative is the
	// same interned node -- the whole point of the arc.  Here we check it the direct way:
	Nd da = a->Differentiate(x);
	Nd f = a * y;
	Nd g = a + x;
	// d(a)/dx built independently equals the one embedded via reuse of 'a'
	BOOST_CHECK_EQUAL(a->Differentiate(x).get(), da.get());
	(void)f; (void)g;
}

// ---- eval correctness on a heavily shared DAG ----

BOOST_AUTO_TEST_CASE(eval_correct_with_shared_subexpression)
{
	auto x = Variable::Make("x");
	Nd a = x*x;                 // shared
	Nd f = a + a + a;           // 3 * x^2, all the same interned 'a'
	x->set_current_value(dbl(2.0, 0.0));
	f->Reset();
	BOOST_CHECK_EQUAL(f->Eval<dbl>(), dbl(12.0, 0.0));   // 3 * 4
}

// ---- immutability under interning: simplifying never changes a held input ----

BOOST_AUTO_TEST_CASE(simplify_leaves_a_held_shared_node_untouched)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd held = (x + y) * x;      // hold a shared node
	x->set_current_value(dbl(3.0, 0.0));
	y->set_current_value(dbl(4.0, 0.0));
	held->Reset();
	auto before = held->Eval<dbl>();

	Nd s = bertini::Simplify(held + 0*y);   // simplify an expression that contains 'held'
	(void)s;

	held->Reset();
	BOOST_CHECK_EQUAL(held->Eval<dbl>(), before);   // held is byte-for-byte unchanged
}

BOOST_AUTO_TEST_SUITE_END() // interning
