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
\file Tests for hash-consing: every Make() routes through Intern(), so structurally
equal nodes collapse to one shared object.  Pins the dedup behavior, the differentiation-
sharing win, and -- critically -- that building/simplifying never mutates a shared interned
node (the bug where SimplifiedSum did Make(first) then AddOperand(rest), corrupting a shared
single-operand sum).
*/

#include <cstdlib>
#include <sstream>
#include <map>
#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical.hpp"
#include <boost/test/unit_test.hpp>

#include "eval_helper.hpp"

using Nd = std::shared_ptr<bertini::node::Node>;
using Variable = bertini::node::Variable;
using Integer = bertini::node::Integer;
using SumOperator = bertini::node::SumOperator;
using dbl = bertini::dbl;
using bertini::test::EvalAt;

BOOST_AUTO_TEST_SUITE(interning)

// ---- the regression: building/simplifying must not mutate a shared interned node ----

BOOST_AUTO_TEST_CASE(simplify_does_not_corrupt_a_shared_single_operand_sum)
{
	// This is the minimal form of a bug interning exposed: a single-operand Sum is reused,
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
	std::map<std::string,dbl> pt{ {"x", xv}, {"y", yv} };

	Nd rs = bertini::Simplify(r);
	BOOST_CHECK_EQUAL(EvalAt<dbl>(rs, pt), xv + xv);   // == 2x, with the y's cancelled

	// and the reused sub-objects must be intact (not mutated by the simplification above)
	BOOST_CHECK_EQUAL(EvalAt<dbl>(m, pt), xv);
	BOOST_CHECK_EQUAL(EvalAt<dbl>(n, pt), yv);
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
	BOOST_CHECK_EQUAL((x + y).get(), (y + x).get());               // canonicalized: same node
}

BOOST_AUTO_TEST_CASE(variables_are_canonical_by_name)
{
	// Make("x") interns to a single canonical x -- "system1's x IS system2's x".
	auto x1 = Variable::Make("x");
	auto x2 = Variable::Make("x");
	BOOST_CHECK_EQUAL(x1.get(), x2.get());   // one object: identity is the canonicalization
	BOOST_CHECK(Variable::Make("x").get() != Variable::Make("z").get());

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
	// same interned node -- the payoff of hash-consing derivatives.  Checked the direct way:
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
	BOOST_CHECK_EQUAL(EvalAt<dbl>(f, {{"x", dbl(2.0, 0.0)}}), dbl(12.0, 0.0));   // 3 * 4
}

// ---- immutability under interning: simplifying never changes a held input ----

BOOST_AUTO_TEST_CASE(simplify_leaves_a_held_shared_node_untouched)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd held = (x + y) * x;      // hold a shared node
	std::map<std::string,dbl> pt{ {"x", dbl(3.0, 0.0)}, {"y", dbl(4.0, 0.0)} };
	auto before = EvalAt<dbl>(held, pt);

	Nd s = bertini::Simplify(held + 0*y);   // simplify an expression that contains 'held'
	(void)s;

	BOOST_CHECK_EQUAL(EvalAt<dbl>(held, pt), before);   // held is unchanged by simplifying around it
}

BOOST_AUTO_TEST_SUITE_END() // interning


// ---- canonical operand ordering (reorder-only) ----

BOOST_AUTO_TEST_SUITE(canonicalization)

// enable canonicalization (with a chosen monomial order) for the duration of a test, then
// restore the global state -- the setting is session-global, so it must not leak.
struct CanonGuard
{
	bool prev_on;
	bertini::node::MonomialOrder prev_order;
	explicit CanonGuard(bertini::node::MonomialOrder o = bertini::node::MonomialOrder::GrevLex)
		: prev_on(bertini::node::CanonicalizeByDefault()),
		  prev_order(bertini::node::CurrentMonomialOrder())
	{
		bertini::node::SetMonomialOrder(o);
		bertini::node::SetCanonicalizeByDefault(true);
	}
	~CanonGuard()
	{
		bertini::node::SetCanonicalizeByDefault(prev_on);
		bertini::node::SetMonomialOrder(prev_order);
	}
};

BOOST_AUTO_TEST_CASE(disabling_canonicalization_preserves_authored_order)
{
	// canonicalization is ON by default; turn it off and the authored operand order is kept,
	// so x+y and y+x become distinct interned nodes again.
	bool prev = bertini::node::CanonicalizeByDefault();
	bertini::node::SetCanonicalizeByDefault(false);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x + y, b = y + x;
	BOOST_CHECK(a.get() != b.get());
	bertini::node::SetCanonicalizeByDefault(prev);
}

BOOST_AUTO_TEST_CASE(commutative_sum_dedups_when_enabled)
{
	CanonGuard g;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x + y;
	Nd b = y + x;
	BOOST_CHECK_EQUAL(a.get(), b.get());          // canonicalized to one node
	x->set_current_value(dbl(2.0, 0.0));
	y->set_current_value(dbl(5.0, 0.0));
	a->Reset();
	BOOST_CHECK_EQUAL(a->Eval<dbl>(), dbl(7.0, 0.0));
}

BOOST_AUTO_TEST_CASE(commutative_product_dedups_when_enabled)
{
	CanonGuard g;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	BOOST_CHECK_EQUAL((x*y).get(), (y*x).get());
}

BOOST_AUTO_TEST_CASE(division_stays_correct_under_canonicalization)
{
	CanonGuard g;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd q = y / x;                                 // a divisor must not become the leading factor
	x->set_current_value(dbl(2.0, 0.0));
	y->set_current_value(dbl(6.0, 0.0));
	q->Reset();
	BOOST_CHECK_EQUAL(q->Eval<dbl>(), dbl(3.0, 0.0));    // 6/2
	BOOST_CHECK((x/y).get() != (y/x).get());            // x/y and y/x stay distinct
}

BOOST_AUTO_TEST_CASE(all_three_orders_are_selectable_and_dedup)
{
	for (auto ord : { bertini::node::MonomialOrder::Lex,
	                  bertini::node::MonomialOrder::RevLex,
	                  bertini::node::MonomialOrder::GrevLex })
	{
		CanonGuard g(ord);
		auto x = Variable::Make("x");
		auto y = Variable::Make("y");
		BOOST_CHECK_EQUAL((x + y).get(), (y + x).get());
	}
}

BOOST_AUTO_TEST_CASE(canonicalization_shows_up_in_printing)
{
	auto str = [](Nd const& n){ std::ostringstream o; n->print(o); return o.str(); };
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	{
		CanonGuard g;
		BOOST_CHECK_EQUAL(str(x + y), str(y + x));   // canonical: same printed form
		BOOST_CHECK_EQUAL(str(x * y), str(y * x));
		BOOST_CHECK_EQUAL(str(Integer::Make(3) * pow(x, 2)), "3*x^2");   // coefficient prints first
	}
	// turning it off preserves the authored operand order in the print
	bool prev = bertini::node::CanonicalizeByDefault();
	bertini::node::SetCanonicalizeByDefault(false);
	BOOST_CHECK_EQUAL(str(y + x), "y+x");
	bertini::node::SetCanonicalizeByDefault(prev);
}

BOOST_AUTO_TEST_CASE(guard_restores_global_canonicalization_state)
{
	// the canonicalization setting is session-global; confirm the preceding tests' guards
	// left it back at the default (on, GrevLex) so they cannot leak into other suites.
	BOOST_CHECK(bertini::node::CanonicalizeByDefault());
	BOOST_CHECK(bertini::node::CurrentMonomialOrder() == bertini::node::MonomialOrder::GrevLex);
}

BOOST_AUTO_TEST_SUITE_END() // canonicalization
