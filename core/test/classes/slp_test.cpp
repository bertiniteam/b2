#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"
#include "bertini2/system/start_systems.hpp"

#include <set>
#include <iostream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <atomic>
#include <gmp.h>

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

	bertini::start_system::TotalDegreeBinomial start(sys);

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

	std::string PrintOfNode(Nd const& n)
	{
		std::ostringstream oss;
		n->print(oss);
		return oss.str();
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
	using bertini::node::MultOperator;
	using bertini::node::IntegerPowerOperator;

	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	std::cout << "\nCSE_TABLE_BEGIN\n";
	std::cout << "| K | polynomial degree (2^K) | SLP slots (fns+Jac) | naive multiplies (2^K-1) |\n";
	std::cout << "|--:|------------------------:|--------------------:|-------------------------:|\n";

	std::size_t prev_slots = 0;
	for (int K = 1; K <= 16; ++K)
	{
		Nd e = x + y;
		for (int i = 0; i < K; ++i)
			e = e * e;          // squaring a repeated base folds: (x+y)^(2^(i+1))

		// The chain no longer builds an exponential binary tree.  Power-folding in
		// CanonicalizeNaryOperands collapses e*e all the way to a single IntegerPower of (x+y)
		// with exponent 2^K -- possibly inside a 1-operand MultOperator wrapper (which is free in
		// the SLP).  Unwrap that wrapper, then assert the folded form.
		Nd inner = e;
		if (auto mo = std::dynamic_pointer_cast<MultOperator const>(inner))
			if (mo->NumOperands() == 1)
				inner = mo->Operands()[0];
		auto ip = std::dynamic_pointer_cast<IntegerPowerOperator const>(inner);
		BOOST_REQUIRE(ip);                                                  // folded to one power
		BOOST_CHECK_EQUAL(ip->exponent(), 1 << K);                         // exponent 2^K
		BOOST_CHECK_EQUAL(PrintOfNode(ip->Operand()), std::string("x+y")); // ...of the (x+y) base

		bertini::System sys;
		sys.AddVariableGroup(bertini::VariableGroup{x, y});
		sys.AddFunction(e);
		const auto slots = SlpSlots(sys);

		const std::size_t naive = (std::size_t{1} << K) - 1;   // multiplies a tree-walk would emit
		std::cout << "| " << K << " | " << (std::size_t{1} << K) << " | "
		          << slots << " | " << naive << " |\n";

		// The compiled program is O(K), not O(2^K): exponentiation-by-squaring emits ~a constant
		// number of instructions per level, so slots grow linearly and stay far below the naive
		// degree-many multiplies.  This is the fold + SLP win the earlier bisection motivated.
		if (K >= 5)   // slots are linear (+~5/level) while naive multiplies double; they cross at K=5
			BOOST_CHECK_LT(slots, naive);
		if (K > 1)
			BOOST_CHECK_LT(slots - prev_slots, std::size_t{16});   // bounded growth per level
		prev_slots = slots;
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

// Regression: the tier dispatch builds complex temporaries when promoting a real operand to complex
// (Power/sqrt/log/...).  Boost inits a fresh mpc at the thread-default precision, which on the eval
// path is not guaranteed to be the working precision -- if it is 0, mpc_init2 aborts.  Eval pins the
// thread precision to the working precision first.  Exercises the promotion paths at mp precision:
// x/y's Jacobian (-x/y^2) promotes a real exponent constant, and sqrt/log/x^3 are the escape ops.
BOOST_AUTO_TEST_CASE(mp_promotion_paths_evaluate_correctly)
{
	bertini::DefaultPrecision(40);
	auto slp = SLP(ParseSLPSys("function f; variable_group x, y; f = x/y;"));
	slp.precision(40);

	Vec<complex_mp> pt(2);
	pt(0) = complex_mp(6); pt(1) = complex_mp(2);
	slp.Eval(pt);
	auto f = slp.GetFuncVals<complex_mp>();
	auto J = slp.GetJacobian<complex_mp>();         // -x/y^2 promotes a real exponent constant

	BOOST_CHECK(abs(f(0)   - complex_mp(3)) < 1e-30);             // 6/2 = 3
	BOOST_CHECK(abs(J(0,0) - complex_mp(1) / complex_mp(2)) < 1e-30); // d(x/y)/dx = 1/y = 1/2

	auto g = SLP(ParseSLPSys("function f; variable_group x; f = sqrt(x) + log(x) + x^3;"));
	g.precision(40);
	Vec<complex_mp> q(1); q(0) = complex_mp(4);
	g.Eval(q);
	auto gf = g.GetFuncVals<complex_mp>();
	using std::sqrt; using std::log;
	BOOST_CHECK(abs(gf(0) - (sqrt(complex_mp(4)) + log(complex_mp(4)) + complex_mp(64))) < 1e-30);
}

// Opt-in A/B speedup benchmark (ADR-0034 gate).  Skipped unless BERTINI_SLP_BENCH is set, so it adds
// no time to normal runs.  Builds the SAME real-coefficient-heavy system twice -- once with tiers
// forced off (all-complex baseline) and once on -- and times N evaluations of function + Jacobian at
// high mpfr precision, where mpfr-complex arithmetic dominates.  Run with:
//   BERTINI_SLP_BENCH=1 ./build/core/test_classes --run_test=SLP_tiered_numtype/tier_speedup_benchmark
BOOST_AUTO_TEST_CASE(tier_speedup_benchmark)
{
	if (!std::getenv("BERTINI_SLP_BENCH")) { BOOST_CHECK(true); return; }

	using real_mp = bertini::real_mp;
	// Generate an nv-variable, nv-function dense integer-coefficient system (degree up to `deg`),
	// so the benchmark can scale: bigger systems do more arithmetic per eval, shrinking the fixed
	// per-eval overhead (allocations, output extraction) relative to the tiered arithmetic.
	const char* nvenv = std::getenv("BERTINI_SLP_BENCH_NVARS");
	const int nv  = nvenv ? std::atoi(nvenv) : 4;
	const int deg = 5;
	// zero-padded names (x00, x01, ...) so no name is a prefix of another (x1 vs x10 confuses the parser)
	auto vn = [](int i){ char b[8]; std::snprintf(b, sizeof b, "x%03d", i); return std::string(b); };
	auto fn = [](int i){ char b[8]; std::snprintf(b, sizeof b, "f%03d", i); return std::string(b); };
	auto make_system = [&]() -> std::string {
		std::string s = "function ";
		for (int i = 0; i < nv; ++i) s += fn(i) + (i + 1 < nv ? "," : "; ");
		s += "variable_group ";
		for (int i = 0; i < nv; ++i) s += vn(i) + (i + 1 < nv ? "," : "; ");
		int c = 2;
		for (int i = 0; i < nv; ++i) {
			s += fn(i) + " = ";
			for (int j = 0; j < nv; ++j)  // one power term per variable
				s += std::to_string(c++) + "*" + vn(j) + "^" + std::to_string((j % deg) + 2) + " + ";
			for (int j = 0; j + 1 < nv; ++j)  // cross terms
				s += std::to_string(c++) + "*" + vn(j) + "*" + vn(j + 1) + " + ";
			s += std::to_string(c++) + "; ";
		}
		return s;
	};
	const std::string dense = make_system();

	const char* penv = std::getenv("BERTINI_SLP_BENCH_PREC");
	const char* nenv = std::getenv("BERTINI_SLP_BENCH_N");
	const unsigned prec = penv ? static_cast<unsigned>(std::atoi(penv)) : 256;
	const int N = nenv ? std::atoi(nenv) : 400;

	auto run = [&](bool tiers) -> double {
		bertini::SLPProgram::tiers_enabled_ = tiers;
		bertini::DefaultPrecision(prec);
		auto sys = ParseSLPSys(dense);
		SLP slp(sys);
		slp.precision(prec);

		Vec<complex_mp> pt(nv);
		for (int j = 0; j < nv; ++j) pt(j) = complex_mp(real_mp(j + 2), real_mp(j + 1));
		slp.Eval(pt); (void)slp.GetFuncVals<complex_mp>();  // warm the frozen prologue

		complex_mp sink(0);
		auto t0 = std::chrono::steady_clock::now();
		for (int i = 0; i < N; ++i) {
			pt(0) = complex_mp(real_mp(i % 7 + 2), real_mp(i % 5 + 1));
			slp.Eval(pt);
			auto f = slp.GetFuncVals<complex_mp>();
			auto J = slp.GetJacobian<complex_mp>();
			sink += f(0) + J(0,0);
		}
		auto t1 = std::chrono::steady_clock::now();
		std::cout << "  [tiers=" << tiers << " real_slots=" << slp.NumRealSlots()
		          << " sink=" << sink.real() << "]\n";
		return std::chrono::duration<double>(t1 - t0).count();
	};

	// DOUBLE precision -- where the AMP tracker spends most of its time, and the worst case for
	// per-op dispatch overhead (cheap hardware arithmetic, so the slot_numtype_ branches dominate).
	const char* ndenv = std::getenv("BERTINI_SLP_BENCH_NDBL");
	const int Nd = ndenv ? std::atoi(ndenv) : 100000;
	auto run_dbl = [&](bool tiers) -> double {
		bertini::SLPProgram::tiers_enabled_ = tiers;
		auto sys = ParseSLPSys(dense);
		SLP slp(sys);
		Vec<complex_dbl> pt(nv);
		for (int j = 0; j < nv; ++j) pt(j) = complex_dbl(j + 2, j + 1);
		slp.Eval(pt); (void)slp.GetFuncVals<complex_dbl>();

		complex_dbl sink(0);
		auto t0 = std::chrono::steady_clock::now();
		for (int i = 0; i < Nd; ++i) {
			pt(0) = complex_dbl(i % 7 + 2, i % 5 + 1);
			slp.Eval(pt);
			auto f = slp.GetFuncVals<complex_dbl>();
			auto J = slp.GetJacobian<complex_dbl>();
			sink += f(0) + J(0,0);
		}
		auto t1 = std::chrono::steady_clock::now();
		std::cout << "  [dbl tiers=" << tiers << " real_slots=" << slp.NumRealSlots()
		          << " sink=" << sink.real() << "]\n";
		return std::chrono::duration<double>(t1 - t0).count();
	};
	const double d_off = run_dbl(false);
	const double d_on  = run_dbl(true);

	const double off = run(false);
	const double on  = run(true);
	bertini::SLPProgram::tiers_enabled_ = true;  // restore default

	std::cout << "\n=== SLP tier benchmark (4 fns / 4 vars) ===\n"
	          << "  DOUBLE  (N=" << Nd << " evals of f+J):\n"
	          << "    tiers OFF: " << d_off << " s    tiers ON: " << d_on
	          << " s    speedup: " << (d_off / d_on) << "x\n"
	          << "  MPFR prec=" << prec << " digits (N=" << N << " evals of f+J):\n"
	          << "    tiers OFF: " << off << " s    tiers ON: " << on
	          << " s    speedup: " << (off / on) << "x\n\n";
	BOOST_CHECK(true);
}

// Probe the user's hypothesis: is real_mp * complex_mp actually cheaper than complex*complex, or
// does Boost.Multiprecision promote the real operand to complex first (4 muls, no saving)?  And how
// much of mpfr cost is the multiply vs. number management?  Opt-in via BERTINI_SLP_BENCH.
BOOST_AUTO_TEST_CASE(mpfr_raw_op_microbench)
{
	if (!std::getenv("BERTINI_SLP_BENCH")) { BOOST_CHECK(true); return; }
	using real_mp = bertini::real_mp;
	const char* penv = std::getenv("BERTINI_SLP_BENCH_PREC");
	const unsigned prec = penv ? static_cast<unsigned>(std::atoi(penv)) : 256;
	const char* menv = std::getenv("BERTINI_SLP_BENCH_M");
	const long M = menv ? std::atol(menv) : 1000000;
	bertini::DefaultPrecision(prec);

	complex_mp ac(real_mp("1.2345"), real_mp("5.6789"));
	complex_mp bc(real_mp("9.8765"), real_mp("4.3210"));
	complex_mp zc(0); bertini::Precision(zc, prec);
	real_mp ar("3.14159"), br("2.71828"), zr(0);

	complex_mp sink(0); bertini::Precision(sink, prec);
	real_mp rsink(0);
	auto timeit = [&](auto&& f){
		auto t0 = std::chrono::steady_clock::now();
		for (long i = 0; i < M; ++i) f();
		return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	};

	const double t_cc = timeit([&]{ zc = ac * bc; sink += zc; });  // complex * complex
	const double t_rc = timeit([&]{ zc = ar * bc; sink += zc; });  // real * complex (tier hot path)
	const double t_rr = timeit([&]{ zr = ar * br; rsink += zr; }); // real * real

	std::cout << "\n=== raw mpfr multiply microbench (prec=" << prec << ", M=" << M
	          << ", sink=" << sink.real() << "/" << rsink << ") ===\n"
	          << "  complex*complex: " << t_cc << " s\n"
	          << "  real*complex:    " << t_rc << " s   (rc/cc = " << (t_rc / t_cc) << ")\n"
	          << "  real*real:       " << t_rr << " s   (rr/cc = " << (t_rr / t_cc) << ")\n"
	          << "  => " << (t_rc / t_cc < 0.8 ? "real*complex IS cheaper -- tier helps arithmetic"
	                                           : "real*complex ~ complex*complex -- Boost promotes; no arithmetic win")
	          << "\n\n";
	BOOST_CHECK(true);
}

namespace {
	// A counting GMP allocator: delegates to std malloc/realloc/free (so blocks stay interchangeable
	// with the default GMP allocator) while tallying calls.  Used to quantify mpfr churn per eval.
	std::atomic<long> g_alloc{0}, g_realloc{0}, g_free{0}, g_bytes{0};
	void* counting_alloc(size_t n)                       { ++g_alloc;  g_bytes += static_cast<long>(n); return std::malloc(n); }
	void* counting_realloc(void* p, size_t /*o*/, size_t n){ ++g_realloc; return std::realloc(p, n); }
	void  counting_free(void* p, size_t /*n*/)           { ++g_free;   std::free(p); }
}

// Quantify the mpfr allocation churn of the live eval segment (the user's hypothesis: alloc/dealloc,
// not op-count, is the mp performance limit).  Counts malloc/realloc/free per eval of f+J at mp
// precision.  Opt-in via BERTINI_SLP_BENCH.
BOOST_AUTO_TEST_CASE(mpfr_alloc_churn_per_eval)
{
	if (!std::getenv("BERTINI_SLP_BENCH")) { BOOST_CHECK(true); return; }
	using real_mp = bertini::real_mp;
	const char* penv = std::getenv("BERTINI_SLP_BENCH_PREC");
	const unsigned prec = penv ? static_cast<unsigned>(std::atoi(penv)) : 256;
	bertini::DefaultPrecision(prec);

	const char* sysenv = std::getenv("BERTINI_SLP_BENCH_SYS");
	auto sys_for_slp = ParseSLPSys(sysenv ? sysenv :
		"function f0,f1,f2,f3; variable_group x0,x1,x2,x3; "
		"f0 = 3*x0^4 + 5*x1^3 + 7*x2^2 + 11*x3^5 + 2*x0*x1*x2 - 13*x2*x3 + 17*x0 - 19; "
		"f1 = 23*x1^4 + 29*x2^3 + 31*x0^2 + 37*x3 + 41*x0*x2*x3 - 43*x0*x1 + 47; "
		"f2 = 53*x2^4 + 59*x3^3 + 61*x0^2 + 67*x1^5 + 71*x1*x2*x3 - 73*x0*x3 + 79; "
		"f3 = 83*x3^4 + 89*x0^3 + 97*x1^2 + 101*x2^4 + 103*x0*x1*x3 - 107*x1*x2 + 109;");
	auto slp = SLP(sys_for_slp);
	if (std::getenv("BERTINI_SLP_DUMP")) std::cout << "\n----- SLP -----\n" << slp << "\n---------------\n";
	slp.precision(prec);
	const int nvars = static_cast<int>(slp.NumVariables());
	Vec<complex_mp> pt(nvars);
	for (int j = 0; j < nvars; ++j) pt(j) = complex_mp(real_mp(j + 2), real_mp(j + 1));
	slp.Eval(pt); (void)slp.GetFuncVals<complex_mp>(); (void)slp.GetJacobian<complex_mp>();  // warm prologue

	// install the counting allocator (blocks stay malloc/free-compatible, so mixing is safe)
	void* (*oa)(size_t); void* (*ora)(void*, size_t, size_t); void (*ofr)(void*, size_t);
	mp_get_memory_functions(&oa, &ora, &ofr);
	mp_set_memory_functions(counting_alloc, counting_realloc, counting_free);

	const int N = 200;
	// Pre-build the varying inputs OUTSIDE the counted region -- constructing complex_mp/real_mp
	// allocates, and we must not attribute that to eval.  Assigning a prebuilt value is alloc-free.
	std::vector<complex_mp> inputs;
	for (int i = 0; i < N; ++i) inputs.push_back(complex_mp(real_mp(i % 7 + 2), real_mp(i % 5 + 1)));

	// (a) eval live segment only
	g_alloc = g_realloc = g_free = 0;
	for (int i = 0; i < N; ++i) { pt(0) = inputs[i]; slp.Eval(pt); }
	const double e_malloc = double(g_alloc)/N, e_realloc = double(g_realloc)/N, e_free = double(g_free)/N;
	// (b) eval + output extraction (GetFuncVals + GetJacobian build Eigen Vec/Mat<mpc>)
	g_alloc = g_realloc = g_free = 0;
	for (int i = 0; i < N; ++i) {
		pt(0) = inputs[i];
		slp.Eval(pt); (void)slp.GetFuncVals<complex_mp>(); (void)slp.GetJacobian<complex_mp>();
	}
	const double t_malloc = double(g_alloc)/N, t_realloc = double(g_realloc)/N, t_free = double(g_free)/N;
	mp_set_memory_functions(oa, ora, ofr);         // restore

	std::cout << "\n=== mpfr alloc churn (prec=" << prec << ", N=" << N << ", 4 fns/4 vars) ===\n"
	          << "  eval only        : malloc=" << e_malloc << " realloc=" << e_realloc << " free=" << e_free << "\n"
	          << "  eval + extraction: malloc=" << t_malloc << " realloc=" << t_realloc << " free=" << t_free << "\n"
	          << "  extraction adds  : malloc=" << (t_malloc-e_malloc) << " realloc=" << (t_realloc-e_realloc) << "\n\n";
	BOOST_CHECK(true);
}

// How many heap ops does a single raw mpc/mpfr op cost?  Tells us how much of the eval churn is
// Boost's per-op INTERNAL scratch (which only a pooling allocator can fix) vs our own temporaries.
BOOST_AUTO_TEST_CASE(mpfr_alloc_per_raw_op)
{
	if (!std::getenv("BERTINI_SLP_BENCH")) { BOOST_CHECK(true); return; }
	using real_mp = bertini::real_mp;
	const char* penv = std::getenv("BERTINI_SLP_BENCH_PREC");
	const unsigned prec = penv ? static_cast<unsigned>(std::atoi(penv)) : 256;
	bertini::DefaultPrecision(prec);
	const long M = 5000;

	complex_mp a(real_mp("1.5"), real_mp("2.5")), b(real_mp("3.5"), real_mp("4.5")), c(0);
	bertini::Precision(c, prec);

	void* (*oa)(size_t); void* (*ora)(void*, size_t, size_t); void (*ofr)(void*, size_t);
	mp_get_memory_functions(&oa, &ora, &ofr);
	auto measure = [&](const char* label, auto&& op){
		mp_set_memory_functions(counting_alloc, counting_realloc, counting_free);
		g_alloc = g_realloc = g_free = 0;
		for (long i = 0; i < M; ++i) op();
		mp_set_memory_functions(oa, ora, ofr);
		std::cout << "  " << label << ": malloc=" << double(g_alloc)/M
		          << " realloc=" << double(g_realloc)/M << " free=" << double(g_free)/M << "\n";
	};
	std::cout << "\n=== heap ops per raw op (prec=" << prec << ", into a preallocated dest) ===\n";
	real_mp r("2.5"); bertini::Precision(r, prec);
	measure("c = a * b   (complex*complex)    ", [&]{ c = a * b; });
	measure("c = a * a   (SQUARING, aliased)  ", [&]{ c = a * a; });
	measure("c = a + b   (complex+complex)    ", [&]{ c = a + b; });
	measure("c = r * a   (real*complex MIXED) ", [&]{ c = r * a; });
	measure("c = r + a   (real+complex MIXED) ", [&]{ c = r + a; });
	measure("c = a * r   (complex*real MIXED) ", [&]{ c = a * r; });
	measure("c = a       (copy / Assign)      ", [&]{ c = a; });
	measure("c = mul-lambda(a,b) (eval binop) ", [&]{ c = [](auto const& x, auto const& y){ return x*y; }(a, b); });
	measure("c = (-x)-lambda(a)  (eval Negate)", [&]{ c = [](auto const& x){ return -x; }(a); });
	measure("c = pow(a,4)(IntPower)           ", [&]{ c = pow(a, 4); });
	// in-place exponentiation by squaring with one reused scratch -- the proposed IntPower replacement
	measure("c = a^4 by squaring (1 scratch)  ", [&]{
		static thread_local complex_mp base; bertini::Precision(base, prec);
		base = a; c = a; for (int k = 1; k < 4; ++k) c *= base;   // simple repeated multiply (alloc-free?)
	});
	BOOST_CHECK(true);
}

// Is the lowered O(n) repeated-multiply actually faster than the transcendental pow(complex,complex)?
// Find the crossover exponent above which we should stop lowering and keep a general pow.  Opt-in.
BOOST_AUTO_TEST_CASE(power_method_crossover)
{
	if (!std::getenv("BERTINI_SLP_BENCH")) { BOOST_CHECK(true); return; }
	using real_mp = bertini::real_mp;
	const char* penv = std::getenv("BERTINI_SLP_BENCH_PREC");
	const unsigned prec = penv ? static_cast<unsigned>(std::atoi(penv)) : 256;
	bertini::DefaultPrecision(prec);
	complex_mp base(real_mp("1.3"), real_mp("0.7")), basec(0), acc(0), tmp(0), result(0);
	for (auto* p : {&base, &basec, &acc, &tmp, &result}) bertini::Precision(*p, prec);
	basec = base;
	const long M = 5000;
	auto timeit = [&](auto&& f){ auto t0 = std::chrono::steady_clock::now();
		for (long i = 0; i < M; ++i) f();
		return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / M * 1e6; };  // us/op

	std::cout << "\n=== power method crossover (prec=" << prec << " digits, us per op) ===\n";
	std::cout << "   n | repeated-mul | pow(c,int) | pow(c,complex) | mul faster than pow(c,c)?\n";
	for (int n : {2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96}) {
		complex_mp en(n); bertini::Precision(en, prec);
		const double t_mul = timeit([&]{ acc = base; for (int k = 2; k <= n; ++k) { tmp = acc * basec; acc.swap(tmp); } result.swap(acc); });
		const double t_pi  = timeit([&]{ result = pow(base, n); });
		const double t_pc  = timeit([&]{ result = pow(base, en); });
		std::cout << "  " << (n<10?" ":"") << n << " | " << t_mul << "      | " << t_pi
		          << "    | " << t_pc << "      | " << (t_mul < t_pc ? "yes" : "NO -- pow wins") << "\n";
	}
	BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END() // SLP_tiered_numtype

