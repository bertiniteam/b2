#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

using bertini::MakeVariable;
using SLP = bertini::StraightLineProgram;
template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
using dbl = bertini::dbl;

BOOST_AUTO_TEST_SUITE(SLP_tests)

bertini::System SimpleTestSystem(){
	std::string str = "function f; variable_group x; f = x+1;";

	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	return sys;
}

BOOST_AUTO_TEST_CASE(can_make_from_system)
{
	auto sys = SimpleTestSystem();

	auto slp = SLP(sys);
}


BOOST_AUTO_TEST_CASE(evaluate_simple_system)
{
	auto sys = SimpleTestSystem();

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



BOOST_AUTO_TEST_CASE(evaluate_system2)
{
	std::string str = "function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;";
	bertini::System sys;
	bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);

	auto slp = SLP(sys);

	Vec<dbl> values(2);

	values(0) = dbl(0.5); // x = 0.5
	values(1) = dbl(0.1); // y = 0.1




	slp.Eval(values);
	Vec<dbl> f = slp.GetFuncVals<dbl>();
	bertini::Mat<dbl> J = slp.GetJacobian<dbl>();


	// not returned yet -- point_d parVals, vec_d parDer,  mat_d Jp

	dbl x{values(0)}, y{values(1)};

	BOOST_CHECK_EQUAL(f(0), pow(x,2)+pow(y,2)-1.); // x^2+y^2-1
	BOOST_CHECK_EQUAL(f(1), x-y);


	BOOST_CHECK_EQUAL(J(0,0), 2.*x); // df1/dx = 2x
	BOOST_CHECK_EQUAL(J(0,1), 2.*y); // df1/dy = 2y
	BOOST_CHECK_EQUAL(J(1,0), 1.);   // df2/dx = 1
	BOOST_CHECK_EQUAL(J(1,1), -1.);  // df2/dy = -1 
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

BOOST_AUTO_TEST_SUITE_END()
