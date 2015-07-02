#include <boost/test/unit_test.hpp>



#include "system.hpp"
#include "system_parsing.hpp"


using Var = std::shared_ptr<bertini::Variable>;
using Num = std::shared_ptr<bertini::Number>;


BOOST_AUTO_TEST_SUITE(homogenization)


BOOST_AUTO_TEST_CASE(no_homogenization_needed_x)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x;
	
	BOOST_CHECK(f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK( f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));

}




BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_1)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);

	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
}



BOOST_AUTO_TEST_CASE(homogenization_needed_1_minus_t_x_plus_t_1_minus_x)
{
	std::shared_ptr<bertini::Variable> x = std::make_shared<bertini::Variable>("x");
	std::shared_ptr<bertini::Variable> t = std::make_shared<bertini::Variable>("t");
	auto f1 = (1-t)*x + t*(1-x);

	BOOST_CHECK(!f1->IsHomogeneous());

	Var h = std::make_shared<bertini::Variable>("h");
	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));
}



BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_t)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-t;
	
	BOOST_CHECK(f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

}



BOOST_AUTO_TEST_CASE(no_homogenization_needed_x_minus_y_t)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var t = std::make_shared<bertini::Variable>("t");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x-y*t;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

}






BOOST_AUTO_TEST_CASE(homogenization_needed_sphere)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = pow(x,2) + pow(y,2) + pow(z,2)-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(y));
	BOOST_CHECK(!f1->IsHomogeneous(z));
	BOOST_CHECK(!f1->IsHomogeneous(h));
}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadric)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var z = std::make_shared<bertini::Variable>("z");
	Var h = std::make_shared<bertini::Variable>("h");

	auto f1 = x*y+x*z+y*z-1;
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);
	vars.push_back(y);
	vars.push_back(z);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(y));
	BOOST_CHECK(!f1->IsHomogeneous(z));
	BOOST_CHECK(!f1->IsHomogeneous(h));

}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x + 1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));

}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant_wrt_y)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");
	Var h = std::make_shared<bertini::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	std::vector<Var> vars;
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK( f1->IsHomogeneous(y));
}





BOOST_AUTO_TEST_CASE(nothomogeneous_sin_x)
{
	Var x = std::make_shared<bertini::Variable>("x");

	auto f1 = sin(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_cos_x)
{
	Var x = std::make_shared<bertini::Variable>("x");

	auto f1 = cos(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());
}

BOOST_AUTO_TEST_CASE(nothomogeneous_tan_x)
{
	Var x = std::make_shared<bertini::Variable>("x");

	auto f1 = tan(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_exp_x)
{
	Var x = std::make_shared<bertini::Variable>("x");

	auto f1 = exp(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());
}


BOOST_AUTO_TEST_CASE(nothomogeneous_sqrt_x)
{
	Var x = std::make_shared<bertini::Variable>("x");

	auto f1 = sqrt(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());
}



BOOST_AUTO_TEST_CASE(is_homogeneous_sin_0)
{
	Num n = std::make_shared<bertini::Number>("1");

	auto f1 = sin(n);
	
	BOOST_CHECK(f1->IsHomogeneous());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_cos_1)
{
	Num n = std::make_shared<bertini::Number>("1");

	auto f1 = cos(n);
	
	BOOST_CHECK(f1->IsHomogeneous());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_sin_1_plus_1)
{
	Num n = std::make_shared<bertini::Number>("1");

	auto f1 = sin(n + n);
	
	BOOST_CHECK(f1->IsHomogeneous());
}




BOOST_AUTO_TEST_CASE(is_homogeneous_summands_homogeneous)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * 4.12331 * pow(x,1);
	
	auto f1 = a+b;
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsHomogeneous(x));
	BOOST_CHECK(f1->IsHomogeneous(y));

}

BOOST_AUTO_TEST_CASE(not_homogeneous_summands_inhomogeneous)
{
	Var x = std::make_shared<bertini::Variable>("x");
	Var y = std::make_shared<bertini::Variable>("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * 4.12331;
	
	auto f1 = a+b;
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(y));
}


BOOST_AUTO_TEST_SUITE_END()





