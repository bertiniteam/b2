#include <boost/test/unit_test.hpp>

#include "bertini2/system/eval_expression.hpp"
#include "bertini2/function_tree.hpp"

#include <map>
#include <string>

using bertini::EvalExpression;
using bertini::node::Variable;
using bertini::node::Integer;
using Nd = std::shared_ptr<bertini::node::Node>;
using dbl = bertini::dbl;
using mpfr_complex = bertini::mpfr_complex;

BOOST_AUTO_TEST_SUITE(eval_expression_without_a_system)

// Evaluate a bare expression at a point given as a name->value map, with no System owned
// by the caller --- the adapter wraps it in a throwaway one-function System internally.
BOOST_AUTO_TEST_CASE(evaluates_a_bare_polynomial_in_double)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x*x + y; // 2^2 + 5 == 9

	std::map<std::string,dbl> values{ {"x", dbl(2,0)}, {"y", dbl(5,0)} };
	auto result = EvalExpression<dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - dbl(9,0)), 1e-14);
}

// Variables are matched to values by name regardless of the order they appear in the map
// or in the expression --- the adapter discovers and orders variables by name internally.
BOOST_AUTO_TEST_CASE(binds_values_by_name_not_position)
{
	auto a = Variable::Make("a");
	auto b = Variable::Make("b");

	Nd f = a - b; // a=7, b=3 -> 4, even though 'b' is listed first in the map

	std::map<std::string,dbl> values{ {"b", dbl(3,0)}, {"a", dbl(7,0)} };
	auto result = EvalExpression<dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - dbl(4,0)), 1e-14);
}

// A constant expression has no variables; an empty value map is the correct input.
BOOST_AUTO_TEST_CASE(evaluates_a_constant_with_no_variables)
{
	Nd f = Integer::Make(3) * Integer::Make(4); // 12

	std::map<std::string,dbl> values; // empty
	auto result = EvalExpression<dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - dbl(12,0)), 1e-14);
}

// A constant expression evaluates in multiple precision too (the empty-point path).
BOOST_AUTO_TEST_CASE(evaluates_a_constant_in_multiple_precision)
{
	bertini::DefaultPrecision(50);

	Nd f = Integer::Make(3) * Integer::Make(4); // 12

	std::map<std::string,mpfr_complex> values; // empty
	auto result = EvalExpression<mpfr_complex>(f, values);

	BOOST_CHECK(abs(result - mpfr_complex(12)) < 1e-40);
}

// Every variable of the expression must be given a value; a missing one is an error.
BOOST_AUTO_TEST_CASE(missing_variable_value_is_an_error)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x + y;

	std::map<std::string,dbl> values{ {"x", dbl(1,0)} }; // y omitted
	BOOST_CHECK_THROW(EvalExpression<dbl>(f, values), std::runtime_error);
}

// Supplying a value for a name not in the expression is an error (typo guard).
BOOST_AUTO_TEST_CASE(value_for_unknown_variable_is_an_error)
{
	auto x = Variable::Make("x");

	Nd f = x*x;

	std::map<std::string,dbl> values{ {"x", dbl(2,0)}, {"z", dbl(9,0)} }; // z not present
	BOOST_CHECK_THROW(EvalExpression<dbl>(f, values), std::runtime_error);
}

// The same expression evaluates correctly in multiple precision.
BOOST_AUTO_TEST_CASE(evaluates_in_multiple_precision)
{
	bertini::DefaultPrecision(50);

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x*x + y; // 3^2 + 4 == 13

	std::map<std::string,mpfr_complex> values{
		{"x", mpfr_complex(3)}, {"y", mpfr_complex(4)} };
	auto result = EvalExpression<mpfr_complex>(f, values);

	BOOST_CHECK(abs(result - mpfr_complex(13)) < 1e-40);
}

BOOST_AUTO_TEST_SUITE_END()
