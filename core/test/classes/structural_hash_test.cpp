//This file is part of Bertini 2.
//
//structural_hash_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//structural_hash_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with structural_hash_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for Node::Hash() and Node::IsSame() -- the order-sensitive, shallow (operands
compared by pointer) structural predicates that the hash-consing intern table will use.
Nothing is wired into construction yet; these just pin the predicate semantics.
*/

#include <cstdlib>
#include "bertini2/function_tree.hpp"
#include <boost/test/unit_test.hpp>

using Nd = std::shared_ptr<bertini::node::Node>;
using Variable = bertini::node::Variable;
using Integer = bertini::node::Integer;
using Rational = bertini::node::Rational;

BOOST_AUTO_TEST_SUITE(structural_hash)

// ---- value leaves: equal-by-value, distinct-by-value ----

BOOST_AUTO_TEST_CASE(equal_integers_are_same_and_hash_equal)
{
	auto a = Integer::Make(7);
	auto b = Integer::Make(7);
	BOOST_CHECK(a->IsSame(*b));
	BOOST_CHECK(b->IsSame(*a));
	BOOST_CHECK_EQUAL(a->Hash(), b->Hash());
}

BOOST_AUTO_TEST_CASE(distinct_integers_are_not_same)
{
	auto a = Integer::Make(7);
	auto b = Integer::Make(8);
	BOOST_CHECK(!a->IsSame(*b));
	// hashes are allowed to collide in principle, but for these they should differ
	BOOST_CHECK(a->Hash() != b->Hash());
}

BOOST_AUTO_TEST_CASE(integer_and_rational_of_same_value_are_not_same)
{
	auto i = Integer::Make(2);
	auto r = Rational::Make("2", "0");
	BOOST_CHECK(!i->IsSame(*r));   // different dynamic type
}

BOOST_AUTO_TEST_CASE(equal_rationals_are_same)
{
	auto a = Rational::Make("3/4", "0");
	auto b = Rational::Make("3/4", "0");
	BOOST_CHECK(a->IsSame(*b));
	BOOST_CHECK_EQUAL(a->Hash(), b->Hash());
}

// ---- variables are identity (until interned by name in Rung 3) ----

BOOST_AUTO_TEST_CASE(distinct_variables_same_name_are_not_same)
{
	auto x1 = Variable::Make("x");
	auto x2 = Variable::Make("x");
	BOOST_CHECK(!x1->IsSame(*x2));      // different objects, by design (pre-interning)
	BOOST_CHECK(x1->IsSame(*x1));
}

// ---- operators: same structure over shared children -> same ----

BOOST_AUTO_TEST_CASE(sums_over_shared_children_are_same)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x + y;
	Nd b = x + y;       // same x, y objects
	BOOST_CHECK(a->IsSame(*b));
	BOOST_CHECK_EQUAL(a->Hash(), b->Hash());
}

BOOST_AUTO_TEST_CASE(order_sensitive_sum)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x + y;
	Nd b = y + x;       // different operand order -> NOT the same (order-sensitive)
	BOOST_CHECK(!a->IsSame(*b));
}

BOOST_AUTO_TEST_CASE(sum_vs_difference_differ_by_signs)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x + y;
	Nd b = x - y;       // same operands, different signs
	BOOST_CHECK(!a->IsSame(*b));
	BOOST_CHECK(a->Hash() != b->Hash());
}

BOOST_AUTO_TEST_CASE(mult_vs_divide_differ)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = x * y;
	Nd b = x / y;
	BOOST_CHECK(!a->IsSame(*b));
}

BOOST_AUTO_TEST_CASE(shallow_distinct_inner_subtrees_are_not_same)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	// (x+y)*z built twice, but each (x+y) is a fresh object -> outer products are NOT same
	// (shallow / pointer comparison of operands; interning in Rung 3 will collapse these).
	Nd a = (x + y) * z;
	Nd b = (x + y) * z;
	BOOST_CHECK(!a->IsSame(*b));
	// but reusing the *same* inner object makes them same:
	Nd inner = x + y;
	Nd c = inner * z;
	Nd d = inner * z;
	BOOST_CHECK(c->IsSame(*d));
	BOOST_CHECK_EQUAL(c->Hash(), d->Hash());
}

BOOST_AUTO_TEST_CASE(transcendentals_distinguished_by_type)
{
	auto x = Variable::Make("x");
	Nd s = sin(x);
	Nd c = cos(x);
	Nd s2 = sin(x);
	BOOST_CHECK(s->IsSame(*s2));
	BOOST_CHECK_EQUAL(s->Hash(), s2->Hash());
	BOOST_CHECK(!s->IsSame(*c));        // sin vs cos: different concrete unary type
	BOOST_CHECK(s->Hash() != c->Hash());
}

BOOST_AUTO_TEST_CASE(integer_power_folds_exponent)
{
	auto x = Variable::Make("x");
	Nd a = pow(x, 2);
	Nd b = pow(x, 2);
	Nd c = pow(x, 3);
	BOOST_CHECK(a->IsSame(*b));
	BOOST_CHECK_EQUAL(a->Hash(), b->Hash());
	BOOST_CHECK(!a->IsSame(*c));        // different integer exponent
}

// ---- hash is structural (recursive) and precision-independent ----

BOOST_AUTO_TEST_CASE(hash_is_structural_and_ignores_working_precision)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	// two structurally-identical trees over the SAME leaves, but with their inner (x+y)
	// built as separate objects, and set to different working precisions.
	Nd f1 = (x + y) * x + Integer::Make(3) * y;
	Nd f2 = (x + y) * x + Integer::Make(3) * y;
	f1->precision(100);
	f2->precision(50);

	// Hash() is recursive over child hashes, so structurally-equal trees hash equal even
	// though their inner subtrees are distinct objects -- and the working precision does
	// not participate.
	BOOST_CHECK_EQUAL(f1->Hash(), f2->Hash());

	// IsSame() is shallow (operands by pointer): the distinct inner (x+y) objects make these
	// NOT the same.  (A hash collision without IsSame is exactly the case the intern table
	// resolves with IsSame; pre-interning it is expected here.)
	BOOST_CHECK(!f1->IsSame(*f2));
}

BOOST_AUTO_TEST_SUITE_END() // structural_hash
