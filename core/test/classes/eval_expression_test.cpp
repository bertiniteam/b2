#include <boost/test/unit_test.hpp>

#include "bertini2/system/eval_expression.hpp"
#include "bertini2/function_tree.hpp"

#include <map>
#include <string>
#include <vector>
#include <functional>

using bertini::EvalExpression;
using bertini::node::Variable;
using bertini::node::Integer;
using Nd = std::shared_ptr<bertini::node::Node>;
using complex_dbl = bertini::complex_dbl;
using complex_mp = bertini::complex_mp;

BOOST_AUTO_TEST_SUITE(eval_expression_without_a_system)

// Evaluate a bare expression at a point given as a name->value map, with no System owned
// by the caller --- the expression's evaluator is compiled internally and memoized on the node.
BOOST_AUTO_TEST_CASE(evaluates_a_bare_polynomial_in_double)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x*x + y; // 2^2 + 5 == 9

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(2,0)}, {"y", complex_dbl(5,0)} };
	auto result = EvalExpression<complex_dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - complex_dbl(9,0)), 1e-14);
}

// Variables are matched to values by name regardless of the order they appear in the map
// or in the expression --- variables are discovered and ordered by name internally.
BOOST_AUTO_TEST_CASE(binds_values_by_name_not_position)
{
	auto a = Variable::Make("a");
	auto b = Variable::Make("b");

	Nd f = a - b; // a=7, b=3 -> 4, even though 'b' is listed first in the map

	std::map<std::string,complex_dbl> values{ {"b", complex_dbl(3,0)}, {"a", complex_dbl(7,0)} };
	auto result = EvalExpression<complex_dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - complex_dbl(4,0)), 1e-14);
}

// A constant expression has no variables; an empty value map is the correct input.
BOOST_AUTO_TEST_CASE(evaluates_a_constant_with_no_variables)
{
	Nd f = Integer::Make(3) * Integer::Make(4); // 12

	std::map<std::string,complex_dbl> values; // empty
	auto result = EvalExpression<complex_dbl>(f, values);

	BOOST_CHECK_SMALL(std::abs(result - complex_dbl(12,0)), 1e-14);
}

// A constant expression evaluates in multiple precision too (the empty-point path).
BOOST_AUTO_TEST_CASE(evaluates_a_constant_in_multiple_precision)
{
	bertini::DefaultPrecision(50);

	Nd f = Integer::Make(3) * Integer::Make(4); // 12

	std::map<std::string,complex_mp> values; // empty
	auto result = EvalExpression<complex_mp>(f, values);

	BOOST_CHECK(abs(result - complex_mp(12)) < 1e-40);
}

// Every variable of the expression must be given a value; a missing one is an error.
BOOST_AUTO_TEST_CASE(missing_variable_value_is_an_error)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x + y;

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(1,0)} }; // y omitted
	BOOST_CHECK_THROW(EvalExpression<complex_dbl>(f, values), std::runtime_error);
}

// Supplying a value for a name not in the expression is an error (typo guard).
BOOST_AUTO_TEST_CASE(value_for_unknown_variable_is_an_error)
{
	auto x = Variable::Make("x");

	Nd f = x*x;

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(2,0)}, {"z", complex_dbl(9,0)} }; // z not present
	BOOST_CHECK_THROW(EvalExpression<complex_dbl>(f, values), std::runtime_error);
}

// With strict=false, a value supplied for a name not in the expression is ignored: the
// expression is simply constant with respect to it (so one superset point serves an
// expression and its variable-dropping derivatives).
BOOST_AUTO_TEST_CASE(strict_false_ignores_unknown_variable)
{
	auto x = Variable::Make("x");

	Nd f = x*x; // x=3 -> 9; the stray 'z' is ignored

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(3,0)}, {"z", complex_dbl(7,0)} };
	auto result = EvalExpression<complex_dbl>(f, values, /*strict=*/false);

	BOOST_CHECK_SMALL(std::abs(result - complex_dbl(9,0)), 1e-14);
}

// strict=false relaxes only the extra-name guard: a variable the expression DOES depend on
// still must be supplied a value, or evaluation throws.
BOOST_AUTO_TEST_CASE(strict_false_still_requires_needed_variable)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x + y;

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(3,0)} }; // y omitted
	BOOST_CHECK_THROW(EvalExpression<complex_dbl>(f, values, /*strict=*/false), std::runtime_error);
}

// The same expression evaluates correctly in multiple precision.
BOOST_AUTO_TEST_CASE(evaluates_in_multiple_precision)
{
	bertini::DefaultPrecision(50);

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	Nd f = x*x + y; // 3^2 + 4 == 13

	std::map<std::string,complex_mp> values{
		{"x", complex_mp(3)}, {"y", complex_mp(4)} };
	auto result = EvalExpression<complex_mp>(f, values);

	BOOST_CHECK(abs(result - complex_mp(13)) < 1e-40);
}

// The compiled evaluator is memoized on the (immutable, hash-consed) expression node: the first
// evaluation compiles and caches it; subsequent evaluations reuse the same cached program rather
// than recompiling a fresh one each call.
BOOST_AUTO_TEST_CASE(memoizes_the_compiled_program_on_the_node)
{
	auto x = Variable::Make("x");
	Nd f = x*x + Integer::Make(1);

	BOOST_CHECK(f->EvalProgram() == nullptr);            // nothing compiled yet

	std::map<std::string,complex_dbl> values{ {"x", complex_dbl(3,0)} };
	auto r1 = EvalExpression<complex_dbl>(f, values);            // x=3 -> 10
	BOOST_CHECK_SMALL(std::abs(r1 - complex_dbl(10,0)), 1e-14);

	auto cached = f->EvalProgram();
	BOOST_CHECK(cached != nullptr);                      // first eval compiled + cached the program

	values["x"] = complex_dbl(5,0);
	auto r2 = EvalExpression<complex_dbl>(f, values);            // x=5 -> 26, still correct
	BOOST_CHECK_SMALL(std::abs(r2 - complex_dbl(26,0)), 1e-14);
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
	auto E1 = [&](Nd f, double xv){ return EvalExpression<complex_dbl>(f, {{"x", complex_dbl(xv,0)}}); };
	auto E2 = [&](Nd f, double xv, double yv){ return EvalExpression<complex_dbl>(f, {{"x", complex_dbl(xv,0)}, {"y", complex_dbl(yv,0)}}); };

	BOOST_CHECK_SMALL(std::abs(E2(x + y + Integer::Make(3), 2, 5) - complex_dbl(10)), tol); // Sum
	BOOST_CHECK_SMALL(std::abs(E2(x - y - Integer::Make(1), 5, 2) - complex_dbl(2)),  tol); // Subtract
	BOOST_CHECK_SMALL(std::abs(E2(x * y * Integer::Make(2), 3, 4) - complex_dbl(24)), tol); // Multiply
	BOOST_CHECK_SMALL(std::abs(E2(x / y, 6, 2) - complex_dbl(3)), tol);                     // Divide
	BOOST_CHECK_SMALL(std::abs(E1(-x, 3) - complex_dbl(-3)), tol);                          // Negate
	BOOST_CHECK_SMALL(std::abs(E1(pow(x, 3), 2) - complex_dbl(8)), tol);                    // IntPower
	BOOST_CHECK_SMALL(std::abs(E2(pow(x, y), 2, 3) - complex_dbl(8)), tol);                 // Power (variable exponent)
	BOOST_CHECK_SMALL(std::abs(E1(sqrt(x), 4) - complex_dbl(2)), tol);                      // Sqrt
	BOOST_CHECK_SMALL(std::abs(E1(exp(x), 0) - complex_dbl(1)), tol);                       // Exp
	BOOST_CHECK_SMALL(std::abs(E1(log(x), 1) - complex_dbl(0)), tol);                       // Log
	BOOST_CHECK_SMALL(std::abs(E1(sin(x), 0) - complex_dbl(0)), tol);                       // Sin
	BOOST_CHECK_SMALL(std::abs(E1(cos(x), 0) - complex_dbl(1)), tol);                       // Cos
	BOOST_CHECK_SMALL(std::abs(E1(tan(x), 0) - complex_dbl(0)), tol);                       // Tan
	BOOST_CHECK_SMALL(std::abs(E1(asin(x), 0) - complex_dbl(0)), tol);                      // Asin
	BOOST_CHECK_SMALL(std::abs(E1(acos(x), 1) - complex_dbl(0)), tol);                      // Acos
	BOOST_CHECK_SMALL(std::abs(E1(atan(x), 0) - complex_dbl(0)), tol);                      // Atan
	BOOST_CHECK_SMALL(std::abs(E1(Pi() * x, 1) - complex_dbl(3.141592653589793, 0)), 1e-12); // Pi
}

// Derivatives produced by Differentiate compile and evaluate correctly through the SLP -- the
// focused successor to the old differentiate-then-node-eval matrix.  The point map is filtered to
// each (simplified) derivative's actual variables.
BOOST_AUTO_TEST_CASE(derivatives_compile_and_evaluate_through_the_slp)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	const double tol = 1e-12;
	auto evald = [](Nd d, std::map<std::string,complex_dbl> known){
		std::map<std::string,complex_dbl> m;
		for (auto const& v : bertini::node::GatherVariables(d)) m[v->name()] = known.at(v->name());
		return EvalExpression<complex_dbl>(d, m);
	};
	std::map<std::string,complex_dbl> pt{ {"x", complex_dbl(2,0)}, {"y", complex_dbl(5,0)} };

	BOOST_CHECK_SMALL(std::abs(evald((x*x)->Differentiate(x), pt)   - complex_dbl(4)),   tol); // d/dx x^2 = 2x, x=2 -> 4
	BOOST_CHECK_SMALL(std::abs(evald((x*y)->Differentiate(x), pt)   - complex_dbl(5)),   tol); // d/dx x*y = y, y=5 -> 5
	BOOST_CHECK_SMALL(std::abs(evald(pow(x,3)->Differentiate(x), pt)- complex_dbl(12)),  tol); // d/dx x^3 = 3x^2, x=2 -> 12
	BOOST_CHECK_SMALL(std::abs(evald(sin(x)->Differentiate(x), pt)  - complex_dbl(std::cos(2.0))), tol); // d/dx sin = cos
	BOOST_CHECK_SMALL(std::abs(evald(exp(x)->Differentiate(x), pt)  - complex_dbl(std::exp(2.0))), tol); // d/dx exp = exp
	BOOST_CHECK_SMALL(std::abs(evald((x/y)->Differentiate(x), pt)   - complex_dbl(0.2)), tol); // d/dx x/y = 1/y, y=5 -> 0.2
}

// Power-folding demonstration + regression guard.  For each expression we compile the SLP with
// the fold OFF (pre-fold canonical form: sort-only) and ON, print BOTH tapes, and assert that
// (a) folding never grows the program and shrinks it when there are repeated factors, and
// (b) the evaluated value is identical with and without the fold.
namespace {
	// Compile f (over the given variables) and return its SLP memory-slot count, printing the tape.
	std::size_t TapeAndSlots(std::shared_ptr<bertini::node::Node> const& f,
	                         bertini::VariableGroup const& vars, char const* label)
	{
		bertini::System sys;
		sys.AddFunction(f);
		sys.AddVariableGroup(vars);
		bertini::StraightLineProgram slp(sys);
		std::cerr << "\n----- " << label << " -----\n";
		std::cerr << "expression prints as: "; f->print(std::cerr); std::cerr << "\n";
		std::cerr << slp << "memory slots (fn+Jac) = " << slp.NumMemorySlots() << "\n";
		return slp.NumMemorySlots();
	}

	complex_dbl EvalAtFixedPoint(std::shared_ptr<bertini::node::Node> const& f)
	{
		std::map<std::string,complex_dbl> all{
			{"x", complex_dbl(0.6, -0.2)},
			{"y", complex_dbl(-1.1, 0.4)},
			{"z", complex_dbl(0.3, 0.7)} };
		std::map<std::string,complex_dbl> pt;    // only the variables this expression actually uses
		for (auto const& v : bertini::node::GatherVariables(f))
			pt[v->name()] = all.at(v->name());
		return EvalExpression<complex_dbl>(f, pt);
	}
}

BOOST_AUTO_TEST_CASE(power_fold_before_after_tape)
{
	using bertini::node::SetPowerFoldByDefault;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	bertini::VariableGroup vars{x, y, z};

	// (name, builder).  The builder runs under each fold setting so the expression is constructed
	// with that canonicalization; leaves (x,y,z) are shared, the products differ.
	struct Case { char const* name; std::function<std::shared_ptr<bertini::node::Node>()> build; };
	std::vector<Case> cases{
		{ "x*x*y",                       [&]{ return x*x*y; } },
		{ "x*x*x*y*y + x*y*y*x",         [&]{ return x*x*x*y*y + x*y*y*x; } },
		{ "y*x*z*x  (reordered square)", [&]{ return y*x*z*x; } },
		{ "(x+y)*(x+y)*(x+y)*z",         [&]{ return (x+y)*(x+y)*(x+y)*z; } },
	};

	for (auto const& c : cases)
	{
		SetPowerFoldByDefault(false);
		auto before = c.build();
		const auto before_slots = TapeAndSlots(before, vars, (std::string("BEFORE fold: ") + c.name).c_str());
		const auto before_val   = EvalAtFixedPoint(before);

		SetPowerFoldByDefault(true);
		auto after = c.build();
		const auto after_slots = TapeAndSlots(after, vars, (std::string("AFTER  fold: ") + c.name).c_str());
		const auto after_val   = EvalAtFixedPoint(after);

		// folding never grows the compiled program ...
		BOOST_CHECK_LE(after_slots, before_slots);
		// ... and never changes the value it computes.
		BOOST_CHECK_SMALL(std::abs(after_val - before_val), 1e-12);
	}

	// Every case here has repeated factors, so the fold strictly shrinks at least one of them.
	SetPowerFoldByDefault(true);   // leave the session in the default state
}

// Regression: SLPProgram::Eval used to read instructions_[ii+3] unconditionally, over-reading the
// instruction tape by one word for a *unary* op (3 words) that is the last instruction -- which is
// the common case, since output wiring ends every program in a trailing Assign. It surfaced as a
// heap-buffer-overflow (ASAN) / SIGSEGV under Guard Malloc for programs whose tape ends exactly on
// an allocation boundary, e.g. the Jacobian of a squared-variable monomial like x*x*y. The value was
// always correct (the stray word is discarded), so only a memory sanitizer/adverse allocator catches
// it -- run this suite under ASAN. Both the bare-expression (EvalExpression) and the normal
// System->StraightLineProgram solve path exercise the same evaluator.
BOOST_AUTO_TEST_CASE(eval_of_trailing_unary_tape_is_memory_safe)
{
	bertini::DefaultPrecision(30);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	std::map<std::string,complex_mp> values{
		{"x", complex_mp(bertini::real_mp("0.6"), bertini::real_mp("-0.2"))},
		{"y", complex_mp(bertini::real_mp("-1.1"), bertini::real_mp("0.4"))} };
	// x^2*y at this point = -0.256 + 0.392i (over-read word never affected the result; check it holds)
	auto v = EvalExpression<complex_mp>(x*x*y, values);
	BOOST_CHECK(abs(v - complex_mp(bertini::real_mp("-0.256"), bertini::real_mp("0.392"))) < 1e-25);
	// exercise several tapes that end in a trailing unary op after a squared-variable factor
	for (Nd f : { Nd(x*x*y), Nd(x*y*y), Nd(x*x/y), Nd(x/(y*y)) })
		EvalExpression<complex_mp>(f, values);
}

BOOST_AUTO_TEST_CASE(system_slp_of_squared_monomial_is_memory_safe)
{
	bertini::DefaultPrecision(30);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	bertini::System sys;
	sys.AddFunction(x*x*y);
	sys.AddVariableGroup(bertini::VariableGroup{x,y});
	bertini::StraightLineProgram slp(sys);      // compiles f + Jacobian; tape ends in a trailing Assign
	bertini::Vec<complex_mp> point(2);
	point << complex_mp(bertini::real_mp("0.6"), bertini::real_mp("-0.2")),
	         complex_mp(bertini::real_mp("-1.1"), bertini::real_mp("0.4"));
	slp.precision(30);
	slp.Eval(point);
	auto fv = slp.GetFuncVals<complex_mp>();
	BOOST_CHECK(abs(fv(0) - complex_mp(bertini::real_mp("-0.256"), bertini::real_mp("0.392"))) < 1e-25);
}

BOOST_AUTO_TEST_SUITE_END()
