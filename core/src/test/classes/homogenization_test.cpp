#include <boost/test/unit_test.hpp>



#include "system.hpp"
#include "system_parsing.hpp"


using Var = std::shared_ptr<bertini::Variable>;

BOOST_AUTO_TEST_SUITE(homogenization)


BOOST_AUTO_TEST_CASE(no_homogenization_needed_x)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x;
	
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
}




BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_y_t)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-1;
	
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
}



BOOST_AUTO_TEST_CASE(homogenization_needed_1_minus_t_x_plus_t_1_minus_x)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");
	auto f1 = (1-t)*x + t*(1-x);


	Var h = std::make_shared<bertini::Variable>("h");
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
}



BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_t)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-t;
	
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
}



BOOST_AUTO_TEST_CASE(no_homogenization_needed_x_minus_y_t)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-y*t;
	
	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
}






BOOST_AUTO_TEST_CASE(homogenization_needed_sphere)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = pow(x,2) + pow(y,2) + pow(z,2)-1;
	
	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadric)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x*y+x*z+y*z-1;
	
	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x + 1;
	
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant_wrt_y)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	std::vector<Var> vars;
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
}





BOOST_AUTO_TEST_SUITE_END()





