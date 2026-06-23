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
// by the caller --- the expression's evaluator is compiled internally and memoized on the node.
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
// or in the expression --- variables are discovered and ordered by name internally.
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

// The compiled evaluator is memoized on the (immutable, hash-consed) expression node: the first
// evaluation compiles and caches it; subsequent evaluations reuse the same cached program rather
// than recompiling a fresh one each call.
BOOST_AUTO_TEST_CASE(memoizes_the_compiled_program_on_the_node)
{
	auto x = Variable::Make("x");
	Nd f = x*x + Integer::Make(1);

	BOOST_CHECK(f->EvalProgram() == nullptr);            // nothing compiled yet

	std::map<std::string,dbl> values{ {"x", dbl(3,0)} };
	auto r1 = EvalExpression<dbl>(f, values);            // x=3 -> 10
	BOOST_CHECK_SMALL(std::abs(r1 - dbl(10,0)), 1e-14);

	auto cached = f->EvalProgram();
	BOOST_CHECK(cached != nullptr);                      // first eval compiled + cached the program

	values["x"] = dbl(5,0);
	auto r2 = EvalExpression<dbl>(f, values);            // x=5 -> 26, still correct
	BOOST_CHECK_SMALL(std::abs(r2 - dbl(26,0)), 1e-14);
	BOOST_CHECK(f->EvalProgram() == cached);             // and reused the very same cached program
}

// Every operator node compiles and evaluates correctly through the straight-line program.  This
// (with the differentiation case below) is the focused successor to the old node-level evaluation
// matrices: node evaluation is gone, so operator correctness is verified through the one engine.
BOOST_AUTO_TEST_CASE(every_operator_compiles_and_evaluates_through_the_slp)
{
	using bertini::node::Pi;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	const double tol = 1e-12;
	auto E1 = [&](Nd f, double xv){ return EvalExpression<dbl>(f, {{"x", dbl(xv,0)}}); };
	auto E2 = [&](Nd f, double xv, double yv){ return EvalExpression<dbl>(f, {{"x", dbl(xv,0)}, {"y", dbl(yv,0)}}); };

	BOOST_CHECK_SMALL(std::abs(E2(x + y + Integer::Make(3), 2, 5) - dbl(10)), tol); // Sum
	BOOST_CHECK_SMALL(std::abs(E2(x - y - Integer::Make(1), 5, 2) - dbl(2)),  tol); // Subtract
	BOOST_CHECK_SMALL(std::abs(E2(x * y * Integer::Make(2), 3, 4) - dbl(24)), tol); // Multiply
	BOOST_CHECK_SMALL(std::abs(E2(x / y, 6, 2) - dbl(3)), tol);                     // Divide
	BOOST_CHECK_SMALL(std::abs(E1(-x, 3) - dbl(-3)), tol);                          // Negate
	BOOST_CHECK_SMALL(std::abs(E1(pow(x, 3), 2) - dbl(8)), tol);                    // IntPower
	BOOST_CHECK_SMALL(std::abs(E2(pow(x, y), 2, 3) - dbl(8)), tol);                 // Power (variable exponent)
	BOOST_CHECK_SMALL(std::abs(E1(sqrt(x), 4) - dbl(2)), tol);                      // Sqrt
	BOOST_CHECK_SMALL(std::abs(E1(exp(x), 0) - dbl(1)), tol);                       // Exp
	BOOST_CHECK_SMALL(std::abs(E1(log(x), 1) - dbl(0)), tol);                       // Log
	BOOST_CHECK_SMALL(std::abs(E1(sin(x), 0) - dbl(0)), tol);                       // Sin
	BOOST_CHECK_SMALL(std::abs(E1(cos(x), 0) - dbl(1)), tol);                       // Cos
	BOOST_CHECK_SMALL(std::abs(E1(tan(x), 0) - dbl(0)), tol);                       // Tan
	BOOST_CHECK_SMALL(std::abs(E1(asin(x), 0) - dbl(0)), tol);                      // Asin
	BOOST_CHECK_SMALL(std::abs(E1(acos(x), 1) - dbl(0)), tol);                      // Acos
	BOOST_CHECK_SMALL(std::abs(E1(atan(x), 0) - dbl(0)), tol);                      // Atan
	BOOST_CHECK_SMALL(std::abs(E1(Pi() * x, 1) - dbl(3.141592653589793, 0)), 1e-12); // Pi
}

// Derivatives produced by Differentiate compile and evaluate correctly through the SLP -- the
// focused successor to the old differentiate-then-node-eval matrix.  The point map is filtered to
// each (simplified) derivative's actual variables.
BOOST_AUTO_TEST_CASE(derivatives_compile_and_evaluate_through_the_slp)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	const double tol = 1e-12;
	auto evald = [](Nd d, std::map<std::string,dbl> known){
		std::map<std::string,dbl> m;
		for (auto const& v : bertini::node::GatherVariables(d)) m[v->name()] = known.at(v->name());
		return EvalExpression<dbl>(d, m);
	};
	std::map<std::string,dbl> pt{ {"x", dbl(2,0)}, {"y", dbl(5,0)} };

	BOOST_CHECK_SMALL(std::abs(evald((x*x)->Differentiate(x), pt)   - dbl(4)),   tol); // d/dx x^2 = 2x, x=2 -> 4
	BOOST_CHECK_SMALL(std::abs(evald((x*y)->Differentiate(x), pt)   - dbl(5)),   tol); // d/dx x*y = y, y=5 -> 5
	BOOST_CHECK_SMALL(std::abs(evald(pow(x,3)->Differentiate(x), pt)- dbl(12)),  tol); // d/dx x^3 = 3x^2, x=2 -> 12
	BOOST_CHECK_SMALL(std::abs(evald(sin(x)->Differentiate(x), pt)  - dbl(std::cos(2.0))), tol); // d/dx sin = cos
	BOOST_CHECK_SMALL(std::abs(evald(exp(x)->Differentiate(x), pt)  - dbl(std::exp(2.0))), tol); // d/dx exp = exp
	BOOST_CHECK_SMALL(std::abs(evald((x/y)->Differentiate(x), pt)   - dbl(0.2)), tol); // d/dx x/y = 1/y, y=5 -> 0.2
}

BOOST_AUTO_TEST_SUITE_END()
