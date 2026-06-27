#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"
#include "bertini2/system/start_systems.hpp"

#include <set>
#include <iostream>

using Variable = bertini::node::Variable;

using bertini::Operation;
using SLP = bertini::StraightLineProgram;
template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;
using complex_dbl = bertini::complex_dbl;

BOOST_AUTO_TEST_SUITE(SLP_tests)


// set up some systems for testing
bertini::System SingleVariableTestSystem(){
	std::string str = "function f; variable_group x; f = x+1;";

	bertini::System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}


bertini::System TwoVariableTestSystem(){
	std::string str = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";

	bertini::System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}



bertini::System ThreeVariableTestSystem(){
	std::string str = "function f, g, h; variable_group x, y, z; f = x+1; g = y-1; h =z/3;";


	bertini::System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}




bertini::System HomotopyTotalDegreeTestSystem(){
	std::string str = "function f, g, h; variable_group x, y, z; f = x+1; g = y-1; h =z/3;";


	bertini::System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	sys.Homogenize();
	sys.AutoPatch();

	
	using Var = std::shared_ptr<Variable>;

	Var t = Variable::Make("t");

	bertini::start_system::TotalDegree start(sys);

	auto homotopy = (1-t)*sys + t*start;

	return homotopy;
}








// begin the actual tests


// super basic -- tests of arities
BOOST_AUTO_TEST_CASE(operation_arities)
{
	BOOST_CHECK(IsUnary(Operation::Negate));
	BOOST_CHECK(IsUnary(Operation::Assign));


	BOOST_CHECK(!IsUnary(Operation::Add));
	BOOST_CHECK(!IsUnary(Operation::Subtract));
	BOOST_CHECK(!IsUnary(Operation::Multiply));
	BOOST_CHECK(!IsUnary(Operation::Divide));
	BOOST_CHECK(!IsUnary(Operation::Power));
}






BOOST_AUTO_TEST_CASE(can_make_from_system)
{
	auto sys = SingleVariableTestSystem();

	auto slp = SLP(sys);
}


BOOST_AUTO_TEST_CASE(has_correct_size)
{
	auto sys = SingleVariableTestSystem();

	auto slp = SLP(sys);

	BOOST_CHECK_EQUAL(slp.NumFunctions(), sys.NumNaturalFunctions());
	BOOST_CHECK_EQUAL(slp.NumVariables(), sys.NumVariables());
}



BOOST_AUTO_TEST_CASE(evaluate_simple_system)
{
	auto sys = SingleVariableTestSystem();

	auto slp = SLP(sys);

	Vec<complex_dbl> values(1);

	values(0) = complex_dbl(2.0);
	

	slp.Eval(values);


	Vec<complex_dbl> f = slp.GetFuncVals<complex_dbl>();
	bertini::Mat<complex_dbl> J = slp.GetJacobian<complex_dbl>();

	// x = 2, and the function is f=x+1
	BOOST_CHECK_EQUAL(f(0), 3.);

	//the system is [f] = [x+1] = [1]

	// so J = matrix of partial derivatives
	//    J = [df/dx] = []

	BOOST_CHECK_EQUAL(J(0,0), 1.);
}


BOOST_AUTO_TEST_CASE(number_variables_system2)
{
	bertini::System sys = TwoVariableTestSystem();
	auto slp = SLP(sys);

	BOOST_CHECK_EQUAL(slp.NumVariables(), sys.NumVariables());

}


BOOST_AUTO_TEST_CASE(evaluate_system2)
{
	bertini::System sys = TwoVariableTestSystem();
	auto slp = SLP(sys);

	Vec<complex_dbl> values(2);

	values(0) = complex_dbl(0.5); // x = 0.5
	values(1) = complex_dbl(0.1); // y = 0.1




	slp.Eval(values);
	Vec<complex_dbl> f = slp.GetFuncVals<complex_dbl>();
	bertini::Mat<complex_dbl> J = slp.GetJacobian<complex_dbl>();


	// not returned yet -- point_d parVals, vec_d parDer,  mat_d Jp

	complex_dbl x{values(0)}, y{values(1)};

	BOOST_CHECK_SMALL(abs(f(0) - (pow(x,2)+pow(y,2)-1.)),1e-15); // x^2+y^2-1
	BOOST_CHECK_SMALL(abs(f(1) - (x-y)),1e-15);


	BOOST_CHECK_SMALL(abs(J(0,0) - (2.*x)),1e-15); // df1/dx = 2x
	BOOST_CHECK_SMALL(abs(J(0,1) - (2.*y)),1e-15); // df1/dy = 2y
	BOOST_CHECK_SMALL(abs(J(1,0) - (1.)),1e-15);   // df2/dx = 1
	BOOST_CHECK_SMALL(abs(J(1,1) - (-1.)),1e-15);  // df2/dy = -1
}



BOOST_AUTO_TEST_CASE(evaluate_system2_inplace)
{
	bertini::System sys = TwoVariableTestSystem();
	auto slp = SLP(sys);

	Vec<complex_dbl> values(2);

	values(0) = complex_dbl(0.5); // x = 0.5
	values(1) = complex_dbl(0.1); // y = 0.1




	slp.Eval(values);


	Vec<complex_dbl> f(slp.NumFunctions());
	slp.GetFuncValsInPlace<complex_dbl>(f);


	bertini::Mat<complex_dbl> J(slp.NumFunctions(), slp.NumVariables());
	slp.GetJacobianInPlace<complex_dbl>(J);


	// not returned yet -- point_d parVals, vec_d parDer,  mat_d Jp

	complex_dbl x{values(0)}, y{values(1)};

	BOOST_CHECK_SMALL(abs(f(0) - (pow(x,2)+pow(y,2)-1.)),1e-15); // x^2+y^2-1
	BOOST_CHECK_SMALL(abs(f(1) - (x-y)),1e-15);


	BOOST_CHECK_SMALL(abs(J(0,0) - (2.*x)),1e-15); // df1/dx = 2x
	BOOST_CHECK_SMALL(abs(J(0,1) - (2.*y)),1e-15); // df1/dy = 2y
	BOOST_CHECK_SMALL(abs(J(1,0) - (1.)),1e-15);   // df2/dx = 1
	BOOST_CHECK_SMALL(abs(J(1,1) - (-1.)),1e-15);  // df2/dy = -1
}



// BOOST_AUTO_TEST_CASE(evaluate){
	// Vec<complex_dbl> values(2);

	// values(0) = complex_dbl(2.0);
	// values(1) = complex_dbl(3.0);

	// Vec<complex_dbl> v = sys.Eval(values);
// 	auto J = sys.Jacobian(values);

// 	complex_dbl x1 = 2;
// 	complex_dbl x2 = 3;

// 	BOOST_CHECK_EQUAL(J(0,0), 2*x1*x2*x2);
// 	BOOST_CHECK_EQUAL(J(0,1), x1*x1*2*x2);
// 	BOOST_CHECK_EQUAL(v(0), 36.0);
// }



BOOST_AUTO_TEST_CASE(has_correct_size_for_three_variable)
{
	auto sys = ThreeVariableTestSystem();

	auto slp = SLP(sys);

	BOOST_CHECK_EQUAL(slp.NumFunctions(), sys.NumNaturalFunctions());
	BOOST_CHECK_EQUAL(slp.NumVariables(), sys.NumVariables());
}



BOOST_AUTO_TEST_CASE(evaluate_three_variable_system)
{
	auto sys = ThreeVariableTestSystem();

	auto slp = SLP(sys);

	Vec<complex_dbl> values(3);

	values(0) = complex_dbl(2.0);
	values(1) = complex_dbl(3.0);
	values(2) = complex_dbl(3.0);


	slp.Eval(values);

	Vec<complex_dbl> f = slp.GetFuncVals<complex_dbl>();
	bertini::Mat<complex_dbl> J = slp.GetJacobian<complex_dbl>();


	BOOST_CHECK_EQUAL(f(0), 3.); // f=x+1, x=2 ==> f=3
	BOOST_CHECK_EQUAL(f(1), 2.); // g = y-1
	BOOST_CHECK_EQUAL(f(2), 1.0); // h =z/3

	//the system is [f] = [x+1] = [1]

	// so J = matrix of partial derivatives
	//    J = [df/dx] = []

	BOOST_CHECK_EQUAL(J(0,0), 1.);
}



BOOST_AUTO_TEST_SUITE_END()


// ---- common-subexpression elimination via hash-consing ----
//
// The SLP compiler keys a node->result-slot map by node pointer and compiles each node only
// once.  Because identical subexpressions are hash-consed to a single shared node, the
// compiled program tracks the *DAG* of distinct subexpressions, not the fully-expanded
// expression tree -- so repeated/shared structure is computed once.

BOOST_AUTO_TEST_SUITE(SLP_cse)

namespace {
	using Node = bertini::node::Node;
	using NaryOperator = bertini::node::NaryOperator;
	using Nd = std::shared_ptr<Node>;

	// nodes in the fully-expanded expression TREE: a shared subnode is counted once per place
	// it appears -- the work a naive, non-CSE tree-walk evaluator would do.
	std::size_t ExpandedTreeNodes(Nd const& n)
	{
		if (auto nary = std::dynamic_pointer_cast<NaryOperator const>(n))
		{
			std::size_t c = 1;
			for (auto const& op : nary->Operands())
				c += ExpandedTreeNodes(op);
			return c;
		}
		return 1;
	}

	// DISTINCT nodes (the DAG): what hash-consing + the SLP compiler actually share.
	void CollectDistinct(Nd const& n, std::set<Node const*>& seen)
	{
		if (!seen.insert(n.get()).second)
			return;
		if (auto nary = std::dynamic_pointer_cast<NaryOperator const>(n))
			for (auto const& op : nary->Operands())
				CollectDistinct(op, seen);
	}
	std::size_t DistinctNodes(Nd const& n)
	{
		std::set<Node const*> seen;
		CollectDistinct(n, seen);
		return seen.size();
	}

	std::size_t SlpSlots(bertini::System const& s)
	{
		bertini::StraightLineProgram slp(s);
		return slp.NumMemorySlots();
	}
}

BOOST_AUTO_TEST_CASE(hash_consing_unifies_independently_built_subexpressions)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	// (x+y) built twice, independently, then added: hash-consing makes them one node, so the
	// DAG has a single (x+y) even though the expanded tree repeats it.
	Nd f = (x + y) + (x + y);
	BOOST_CHECK_EQUAL(DistinctNodes(f), 4u);      // x, y, (x+y), the outer sum
	BOOST_CHECK_EQUAL(ExpandedTreeNodes(f), 7u);  // outer + 2*(sum + x + y)
}

BOOST_AUTO_TEST_CASE(cse_benchmark_squaring_chain)
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	std::cout << "\nCSE_TABLE_BEGIN\n";
	std::cout << "| K | expanded tree nodes | distinct nodes (DAG) | SLP slots (fns+Jac) | reduction (tree/DAG) |\n";
	std::cout << "|--:|--------------------:|---------------------:|--------------------:|---------------------:|\n";

	for (int K = 1; K <= 16; ++K)
	{
		Nd e = x + y;
		for (int i = 0; i < K; ++i)
			e = e * e;          // e*e reuses the same node; the DAG grows by one per level

		const auto expanded = ExpandedTreeNodes(e);
		const auto distinct = DistinctNodes(e);

		bertini::System sys;
		sys.AddVariableGroup(bertini::VariableGroup{x, y});
		sys.AddFunction(e);
		const auto slots = SlpSlots(sys);

		std::cout << "| " << K << " | " << expanded << " | " << distinct
		          << " | " << slots << " | " << (expanded / distinct) << "x |\n";

		// the same function is a linear DAG but an exponential tree: hash-consing collapses it
		BOOST_CHECK_EQUAL(distinct, static_cast<std::size_t>(K + 3));   // x, y, (x+y), e_1..e_K
		BOOST_CHECK_EQUAL(expanded, (std::size_t{1} << (K + 2)) - 1);   // a binary tree
		if (K >= 8)
			BOOST_CHECK_LT(slots, expanded);   // the compiled program stays DAG-sized
	}
	std::cout << "CSE_TABLE_END\n" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END() // SLP_cse


// (the SLP-vs-tree oracle suite lived here; removed when the FunctionTree eval method was
// retired -- the SLP is now the sole system evaluator.)


// ---- freeze-set tape partition (ADR-0027) ----
//
// The compiler stably reorders the instruction tape so every constant-only ("frozen")
// instruction precedes every variable-dependent ("live") one.  A point-only re-evaluation
// then skips the frozen prologue (reusing the constants already in memory), and the whole
// tape runs only when the constants are not yet valid for the working precision.

BOOST_AUTO_TEST_SUITE(SLP_freeze_partition)

using bertini::node::Integer;
using complex_mp = bertini::complex_mp;

namespace {
	// f = x + sin(1): sin(1) is a constant unary operation -> a frozen instruction.
	bertini::System ConstantSubexpressionSystem()
	{
		auto x = Variable::Make("x");
		bertini::System sys;
		sys.AddVariableGroup(bertini::VariableGroup{x});
		sys.AddFunction(x + sin(Integer::Make(1)));
		return sys;
	}
}

BOOST_AUTO_TEST_CASE(constant_subexpression_yields_a_frozen_prologue)
{
	auto slp = SLP(ConstantSubexpressionSystem());
	// sin(1) compiles to at least one frozen instruction, so the live segment does not start at 0.
	BOOST_CHECK_GT(slp.FirstLiveInstructionOffset(), 0u);
}

BOOST_AUTO_TEST_CASE(no_constant_operator_means_no_prologue)
{
	auto x = Variable::Make("x");
	bertini::System sys;
	sys.AddVariableGroup(bertini::VariableGroup{x});
	sys.AddFunction(x * x); // every instruction depends on x -> nothing frozen
	auto slp = SLP(sys);
	BOOST_CHECK_EQUAL(slp.FirstLiveInstructionOffset(), 0u);
}

// Skipping the frozen prologue on a point-only change must not corrupt the result: the
// constant stays in memory and the second eval (at a new point) reuses it.
BOOST_AUTO_TEST_CASE(point_only_change_reuses_constants_correctly)
{
	auto slp = SLP(ConstantSubexpressionSystem());

	const double s1 = std::sin(1.0);

	Vec<complex_dbl> p(1);
	p(0) = complex_dbl(2.0);
	slp.Eval(p);
	BOOST_CHECK_CLOSE(slp.GetFuncVals<complex_dbl>()(0).real(), 2.0 + s1, 1e-10);

	p(0) = complex_dbl(5.0); // point-only change: prologue skipped, sin(1) reused
	slp.Eval(p);
	BOOST_CHECK_CLOSE(slp.GetFuncVals<complex_dbl>()(0).real(), 5.0 + s1, 1e-10);

	p(0) = complex_dbl(2.0); // back again
	slp.Eval(p);
	BOOST_CHECK_CLOSE(slp.GetFuncVals<complex_dbl>()(0).real(), 2.0 + s1, 1e-10);
}

// A precision change must invalidate the frozen prologue so the constant is recomputed at the
// new precision.  The tolerance (1e-40) is tighter than the original precision (30 digits): if
// invalidation were broken and sin(1) stayed at 30-digit accuracy, this would fail.
BOOST_AUTO_TEST_CASE(precision_change_recomputes_constants)
{
	auto sys = ConstantSubexpressionSystem();

	bertini::DefaultPrecision(30);
	sys.precision(30);
	Vec<complex_mp> p30(1);
	p30(0) = complex_mp(2);
	auto f30 = sys.Eval(p30);
	BOOST_CHECK(abs(f30(0) - (complex_mp(2) + sin(complex_mp(1)))) < 1e-25);

	bertini::DefaultPrecision(50);
	sys.precision(50);
	Vec<complex_mp> p50(1);
	p50(0) = complex_mp(2);
	auto f50 = sys.Eval(p50);
	// sin(1) recomputed at 50 digits -> accurate well past 30 digits.
	BOOST_CHECK(abs(f50(0) - (complex_mp(2) + sin(complex_mp(1)))) < 1e-40);
}

BOOST_AUTO_TEST_SUITE_END() // SLP_freeze_partition



// ADR-0034: real-valued subexpressions evaluate in the cheaper real banks (NumType::Real), and the
// result must equal the all-complex evaluation.  The rest of the C++ suite is the broad equivalence
// gate (it solves/tracks integer- and rational-coefficient systems through the SLP); here we confirm
// the inference is actually LIVE (not a silent all-Complex no-op) and that the mixed real-coefficient
// x complex-variable path is numerically correct, at double and multiprecision.
BOOST_AUTO_TEST_SUITE(SLP_tiered_numtype)

using complex_mp = bertini::complex_mp;

namespace {
	bertini::System ParseSLPSys(std::string const& str){
		bertini::System sys;
		[[maybe_unused]] bool ok = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
		return sys;
	}
	// f = 3*x^2 + 2 : integer coefficients 3 and 2 are real, x is complex.
	const std::string kRealCoeffSys = "function f; variable_group x; f = 3*x^2 + 2;";
}

// An integer-coefficient system must infer at least one Real slot -- otherwise the tier machinery
// would be a silent no-op and there would be no speedup.
BOOST_AUTO_TEST_CASE(integer_coeff_system_infers_real_slots)
{
	auto slp = SLP(ParseSLPSys(kRealCoeffSys));
	BOOST_CHECK_GT(slp.NumRealSlots(), size_t(0));
}

// Evaluated at a genuinely complex point, the real-coefficient * complex-variable path must produce
// the correct complex value (and Jacobian) -- i.e. the imaginary part propagates through the mixed
// real*complex arithmetic.  x = 1+i:  3*(1+i)^2 + 2 = 3*(2i) + 2 = 2 + 6i;  df/dx = 6x = 6 + 6i.
BOOST_AUTO_TEST_CASE(real_tier_mixed_eval_correct_double)
{
	auto slp = SLP(ParseSLPSys(kRealCoeffSys));
	BOOST_REQUIRE_GT(slp.NumRealSlots(), size_t(0));

	Vec<complex_dbl> x(1); x(0) = complex_dbl(1.0, 1.0);
	slp.Eval(x);
	auto f = slp.GetFuncVals<complex_dbl>();
	auto J = slp.GetJacobian<complex_dbl>();
	BOOST_CHECK_SMALL(abs(f(0)   - complex_dbl(2.0, 6.0)), 1e-13);
	BOOST_CHECK_SMALL(abs(J(0,0) - complex_dbl(6.0, 6.0)), 1e-13);
}

// Same, at multiprecision: the real_mp bank carries the coefficients at the working precision.
BOOST_AUTO_TEST_CASE(real_tier_mixed_eval_correct_mp)
{
	auto slp = SLP(ParseSLPSys(kRealCoeffSys));
	Vec<complex_mp> x(1); x(0) = complex_mp(1, 1);
	slp.Eval(x);
	auto f = slp.GetFuncVals<complex_mp>();
	BOOST_CHECK(abs(f(0) - complex_mp(2, 6)) < 1e-40);
}

BOOST_AUTO_TEST_SUITE_END() // SLP_tiered_numtype

