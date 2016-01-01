//This file is part of Bertini 2.0.
//
//homogenization_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//homogenization_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with homogenization_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  homogenization_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include <boost/test/unit_test.hpp>



#include "system.hpp"
#include "system_parsing.hpp"


using Var = std::shared_ptr<bertini::node::Variable>;
using Float = std::shared_ptr<bertini::node::Float>;

using VariableGroup = bertini::VariableGroup;

BOOST_AUTO_TEST_SUITE(homogenization)


BOOST_AUTO_TEST_CASE(no_homogenization_needed_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = x;
	
	BOOST_CHECK(f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK( f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));


	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_1)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var t = std::make_shared<bertini::node::Variable>("t");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = x-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);

	BOOST_CHECK(f1->IsHomogeneous());
	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));


	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_1_minus_t_x_plus_t_1_minus_x)
{
	std::shared_ptr<bertini::node::Variable> x = std::make_shared<bertini::node::Variable>("x");
	std::shared_ptr<bertini::node::Variable> t = std::make_shared<bertini::node::Variable>("t");
	auto f1 = (1-t)*x + t*(1-x);

	BOOST_CHECK(!f1->IsHomogeneous());

	Var h = std::make_shared<bertini::node::Variable>("h");
	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_x_minus_t)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var t = std::make_shared<bertini::node::Variable>("t");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = x-t;
	
	BOOST_CHECK(f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(no_homogenization_needed_x_minus_y_t)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");
	Var t = std::make_shared<bertini::node::Variable>("t");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = x-y*t;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));
	BOOST_CHECK(!f1->IsHomogeneous(t));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}






BOOST_AUTO_TEST_CASE(homogenization_needed_sphere)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");
	Var z = std::make_shared<bertini::node::Variable>("z");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = pow(x,2) + pow(y,2) + pow(z,2)-1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
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


	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(homogenization_needed_quadric)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");
	Var z = std::make_shared<bertini::node::Variable>("z");
	Var h = std::make_shared<bertini::node::Variable>("h");

	auto f1 = x*y+x*z+y*z-1;
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
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

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));


	BOOST_CHECK(f1->IsPolynomial());
}	





BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var h = std::make_shared<bertini::node::Variable>("h");
	

	auto f1 = pow(x,2) + x + 1;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 2);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var h = std::make_shared<bertini::node::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(x);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 1);
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK(!f1->IsHomogeneous(h));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(homogenization_needed_quadratic_no_constant_wrt_y)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");
	Var h = std::make_shared<bertini::node::Variable>("h");
	

	auto f1 = pow(x,2) + x;
	
	BOOST_CHECK(!f1->IsHomogeneous());

	VariableGroup vars;
	vars.push_back(y);

	f1->Homogenize(vars,h);
	BOOST_CHECK_EQUAL(f1->Degree(h), 0);
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(h));
	BOOST_CHECK( f1->IsHomogeneous(y));

	vars.push_back(h);
	BOOST_CHECK(f1->IsHomogeneous(vars));

	BOOST_CHECK(f1->IsPolynomial());
}





BOOST_AUTO_TEST_CASE(nothomogeneous_sin_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");

	auto f1 = sin(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_cos_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");

	auto f1 = cos(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(nothomogeneous_tan_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");

	auto f1 = tan(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_exp_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");

	auto f1 = exp(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(nothomogeneous_sqrt_x)
{
	Var x = std::make_shared<bertini::node::Variable>("x");

	auto f1 = sqrt(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(is_homogeneous_sin_0)
{
	Float n = std::make_shared<bertini::node::Float>("1");

	auto f1 = sin(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_cos_1)
{
	Float n = std::make_shared<bertini::node::Float>("1");

	auto f1 = cos(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_sin_1_plus_1)
{
	Float n = std::make_shared<bertini::node::Float>("1");

	auto f1 = sin(n + n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(is_homogeneous_summands_homogeneous)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * 4.12331 * pow(x,1);
	
	auto f1 = a+b;
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsHomogeneous(x));
	BOOST_CHECK(f1->IsHomogeneous(y));


	BOOST_CHECK(f1->IsPolynomial());

}

BOOST_AUTO_TEST_CASE(not_homogeneous_summands_inhomogeneous)
{
	Var x = std::make_shared<bertini::node::Variable>("x");
	Var y = std::make_shared<bertini::node::Variable>("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * 4.12331;
	
	auto f1 = a+b;
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(y));

	BOOST_CHECK(f1->IsPolynomial());
}


BOOST_AUTO_TEST_SUITE_END()





