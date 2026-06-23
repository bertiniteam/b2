#include <boost/test/unit_test.hpp>

#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/find.hpp"
#include "bertini2/system/eval_expression.hpp"

#include <map>
#include <sstream>

using bertini::node::Variable;
using bertini::node::Named;
using bertini::node::NamedExpression;
using bertini::node::Find;
using Nd = std::shared_ptr<bertini::node::Node>;
using dbl = bertini::dbl;
using bertini::EvalExpression;

BOOST_AUTO_TEST_SUITE(named_expression)

// A NamedExpression prints as just its name (the expansion is revealed elsewhere, by Describe).
BOOST_AUTO_TEST_CASE(prints_as_its_name)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = Named(x*x + y*y, "a");

	std::stringstream ss; ss << a;
	BOOST_CHECK_EQUAL(ss.str(), "a");
}

// It evaluates to its wrapped expression.
BOOST_AUTO_TEST_CASE(evaluates_to_its_expression)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = Named(x*x + y*y, "a");
	// a at (3,4) = 9 + 16 = 25
	auto v = EvalExpression<dbl>(a, {{"x", dbl(3,0)}, {"y", dbl(4,0)}});
	BOOST_CHECK_SMALL(std::abs(v - dbl(25,0)), 1e-12);
}

// Hash-consed by (expression, name): same -> one node; different name or bare expr -> different.
BOOST_AUTO_TEST_CASE(hash_consed_by_expression_and_name)
{
	auto x = Variable::Make("x");
	Nd e  = x*x;
	Nd a1 = Named(e, "a");
	Nd a2 = Named(e, "a");
	Nd b  = Named(e, "b");

	BOOST_CHECK(a1.get() == a2.get());  // same expression + name -> the same interned node
	BOOST_CHECK(a1.get() != b.get());   // different name -> a different node
	BOOST_CHECK(a1.get() != e.get());   // distinct from the bare expression it wraps
}

// Usable as a subexpression: a named node referenced multiple times evaluates correctly.
BOOST_AUTO_TEST_CASE(usable_as_a_subexpression)
{
	auto x = Variable::Make("x");
	Nd a = Named(x*x, "a");
	Nd f = a*a + a;   // a^2 + a, with a = x^2; at x=2 -> 16 + 4 = 20
	auto v = EvalExpression<dbl>(f, {{"x", dbl(2,0)}});
	BOOST_CHECK_SMALL(std::abs(v - dbl(20,0)), 1e-12);
}

// Find<NamedExpression> discovers the named subexpressions in a tree, including nested ones,
// sorted by name.
BOOST_AUTO_TEST_CASE(find_discovers_named_expressions)
{
	auto x = Variable::Make("x");
	Nd a = Named(x*x, "a");
	Nd b = Named(a + 1, "b");   // nested: b's entry references a
	Nd f = b*b;

	auto found = Find<NamedExpression>(f);
	BOOST_REQUIRE_EQUAL(found.size(), 2u);   // a and b
	BOOST_CHECK_EQUAL(found[0]->name(), "a"); // sorted by name
	BOOST_CHECK_EQUAL(found[1]->name(), "b");
}

// Find<Variable> is GatherVariables' general form: sorted, deduped, descending through Named.
BOOST_AUTO_TEST_CASE(find_variables_descends_through_named)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	Nd a = Named(x*x, "a");
	Nd f = a + y;   // x lives inside the named expression a; y is bare

	auto vars = Find<Variable>(f);
	BOOST_REQUIRE_EQUAL(vars.size(), 2u);
	BOOST_CHECK_EQUAL(vars[0]->name(), "x");
	BOOST_CHECK_EQUAL(vars[1]->name(), "y");
}

BOOST_AUTO_TEST_SUITE_END()
