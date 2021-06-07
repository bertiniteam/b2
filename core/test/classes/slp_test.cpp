#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

using bertini::MakeVariable;
using SLP = bertini::StraightLineProgram;

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


// BOOST_AUTO_TEST_CASE(evaluate){
// 	Vec<dbl> values(2);

// 	values(0) = dbl(2.0);
// 	values(1) = dbl(3.0);

// 	Vec<dbl> v = sys.Eval(values);
// 	auto J = sys.Jacobian(values);

// 	double x1 = 2;
// 	double x2 = 3;

// 	BOOST_CHECK_EQUAL(J(0,0), 2*x1*x2*x2);
// 	BOOST_CHECK_EQUAL(J(0,1), x1*x1*2*x2);
// 	BOOST_CHECK_EQUAL(v(0), 36.0);
// }

BOOST_AUTO_TEST_SUITE_END()
