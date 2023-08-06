#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"
#include "bertini2/system/start_systems.hpp"

using bertini::MakeVariable;
using bertini::Operation;
using SLP = bertini::StraightLineProgram;
template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
using dbl = bertini::dbl;

BOOST_AUTO_TEST_SUITE(SLP_tests)


// set up some systems for testing
bertini::System SingleVariableTestSystem(){
	std::string str = "function f; variable_group x; f = x+1;";

	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}


bertini::System TwoVariableTestSystem(){
	std::string str = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";

	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}



bertini::System ThreeVariableTestSystem(){
	std::string str = "function f, g, h; variable_group x, y, z; f = x+1; g = y-1; h =z/3;";


	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}




bertini::System HomotopyTotalDegreeTestSystem(){
	std::string str = "function f, g, h; variable_group x, y, z; f = x+1; g = y-1; h =z/3;";


	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	sys.Homogenize();
	sys.AutoPatch();

	using bertini::MakeVariable;
	using Variable = bertini::node::Variable;
	using Var = std::shared_ptr<Variable>;

	Var t = MakeVariable("t");

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
