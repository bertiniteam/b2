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
using dbl = bertini::dbl;

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

	Vec<dbl> values(1);

	values(0) = dbl(2.0);
	

	slp.Eval(values);


	Vec<dbl> f = slp.GetFuncVals<dbl>();
	bertini::Mat<dbl> J = slp.GetJacobian<dbl>();

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

	Vec<dbl> values(2);

	values(0) = dbl(0.5); // x = 0.5
	values(1) = dbl(0.1); // y = 0.1




	slp.Eval(values);
	Vec<dbl> f = slp.GetFuncVals<dbl>();
	bertini::Mat<dbl> J = slp.GetJacobian<dbl>();


	// not returned yet -- point_d parVals, vec_d parDer,  mat_d Jp

	dbl x{values(0)}, y{values(1)};

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

	Vec<dbl> values(2);

	values(0) = dbl(0.5); // x = 0.5
	values(1) = dbl(0.1); // y = 0.1




	slp.Eval(values);


	Vec<dbl> f(slp.NumFunctions());
	slp.GetFuncValsInPlace<dbl>(f);


	bertini::Mat<dbl> J(slp.NumFunctions(), slp.NumVariables());
	slp.GetJacobianInPlace<dbl>(J);


	// not returned yet -- point_d parVals, vec_d parDer,  mat_d Jp

	dbl x{values(0)}, y{values(1)};

	BOOST_CHECK_SMALL(abs(f(0) - (pow(x,2)+pow(y,2)-1.)),1e-15); // x^2+y^2-1
	BOOST_CHECK_SMALL(abs(f(1) - (x-y)),1e-15);


	BOOST_CHECK_SMALL(abs(J(0,0) - (2.*x)),1e-15); // df1/dx = 2x
	BOOST_CHECK_SMALL(abs(J(0,1) - (2.*y)),1e-15); // df1/dy = 2y
	BOOST_CHECK_SMALL(abs(J(1,0) - (1.)),1e-15);   // df2/dx = 1
	BOOST_CHECK_SMALL(abs(J(1,1) - (-1.)),1e-15);  // df2/dy = -1
}



// BOOST_AUTO_TEST_CASE(evaluate){
	// Vec<dbl> values(2);

	// values(0) = dbl(2.0);
	// values(1) = dbl(3.0);

	// Vec<dbl> v = sys.Eval(values);
// 	auto J = sys.Jacobian(values);

// 	dbl x1 = 2;
// 	dbl x2 = 3;

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

	Vec<dbl> values(3);

	values(0) = dbl(2.0);
	values(1) = dbl(3.0);
	values(2) = dbl(3.0);


	slp.Eval(values);

	Vec<dbl> f = slp.GetFuncVals<dbl>();
	bertini::Mat<dbl> J = slp.GetJacobian<dbl>();


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


// ---- oracle: the SLP must agree with the legacy tree-walk evaluator ----
//
// Before retiring EvalMethod::FunctionTree, pin that the compiled SLP produces the same
// functions, Jacobian, and time-derivative as evaluating the function trees directly, at the
// same point (up to floating-point reorder rounding).

BOOST_AUTO_TEST_SUITE(SLP_oracle)

namespace {
	void AssertClose(dbl const& a, dbl const& b, double tol)
	{
		BOOST_CHECK_LE(std::abs(a - b), tol * (1.0 + std::abs(a)));
	}

	// evaluate functions + Jacobian both ways and assert agreement
	void OracleCheck(std::string const& sys_text, Vec<dbl> const& point, double tol)
	{
		bertini::System sys(sys_text);
		sys.Differentiate();

		sys.SetEvalMethod(bertini::EvalMethod::FunctionTree);
		Vec<dbl> f_tree = sys.Eval(point);
		Mat<dbl> J_tree = sys.Jacobian(point);

		sys.SetEvalMethod(bertini::EvalMethod::SLP);
		Vec<dbl> f_slp = sys.Eval(point);
		Mat<dbl> J_slp = sys.Jacobian(point);

		BOOST_REQUIRE_EQUAL(f_tree.size(), f_slp.size());
		for (Eigen::Index ii = 0; ii < f_tree.size(); ++ii)
			AssertClose(f_tree(ii), f_slp(ii), tol);
		BOOST_REQUIRE_EQUAL(J_tree.size(), J_slp.size());
		for (Eigen::Index ii = 0; ii < J_tree.rows(); ++ii)
			for (Eigen::Index jj = 0; jj < J_tree.cols(); ++jj)
				AssertClose(J_tree(ii, jj), J_slp(ii, jj), tol);
	}
}

BOOST_AUTO_TEST_CASE(slp_matches_tree_polynomial)
{
	Vec<dbl> pt(3);
	pt << dbl(0.7, 0.3), dbl(-0.4, 0.6), dbl(0.5, -0.2);
	OracleCheck("variable_group x,y,z; function f1,f2; f1 = x^2*y - z*x; f2 = x*y*z + 1;", pt, 1e-12);
}

BOOST_AUTO_TEST_CASE(slp_matches_tree_transcendental)
{
	Vec<dbl> pt(2);
	pt << dbl(0.6, -0.2), dbl(-0.3, 0.4);
	OracleCheck("variable_group x,y; function f; f = sin(x)*cos(y) + exp(x*y) - x/y;", pt, 1e-12);
}

BOOST_AUTO_TEST_CASE(slp_matches_tree_high_degree_shared_subexpressions)
{
	Vec<dbl> pt(2);
	pt << dbl(0.5, 0.1), dbl(0.4, -0.3);
	OracleCheck("variable_group x,y; function f1,f2; f1 = (x^2+y^2)^3; f2 = (x^2+y^2)*x*y;", pt, 1e-11);
}

BOOST_AUTO_TEST_SUITE_END() // SLP_oracle
