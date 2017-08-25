//This file is part of Bertini 2.
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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

//  homogenization_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include <boost/test/unit_test.hpp>



#include "bertini2/system/system.hpp"



BOOST_AUTO_TEST_SUITE(homogenization)

using mpfr_float = bertini::mpfr_float;
using Var = std::shared_ptr<bertini::node::Variable>;
using Float = std::shared_ptr<bertini::node::Float>;
using VariableGroup = bertini::VariableGroup;

using bertini::MakeVariable;
using bertini::MakeFloat;


BOOST_AUTO_TEST_CASE(no_homogenization_needed_x)
{
	Var x = MakeVariable("x");
	Var h = MakeVariable("h");

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
	Var x = MakeVariable("x");
	Var t = MakeVariable("t");
	Var h = MakeVariable("h");

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
	std::shared_ptr<bertini::node::Variable> x = MakeVariable("x");
	std::shared_ptr<bertini::node::Variable> t = MakeVariable("t");
	auto f1 = (1-t)*x + t*(1-x);

	BOOST_CHECK(!f1->IsHomogeneous());

	Var h = MakeVariable("h");
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
	Var x = MakeVariable("x");
	Var t = MakeVariable("t");
	Var h = MakeVariable("h");

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
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");
	Var t = MakeVariable("t");
	Var h = MakeVariable("h");

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
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");
	Var z = MakeVariable("z");
	Var h = MakeVariable("h");

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
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");
	Var z = MakeVariable("z");
	Var h = MakeVariable("h");

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
	Var x = MakeVariable("x");
	Var h = MakeVariable("h");
	

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
	Var x = MakeVariable("x");
	Var h = MakeVariable("h");
	

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
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");
	Var h = MakeVariable("h");
	

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
	Var x = MakeVariable("x");

	auto f1 = sin(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_cos_x)
{
	Var x = MakeVariable("x");

	auto f1 = cos(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(nothomogeneous_tan_x)
{
	Var x = MakeVariable("x");

	auto f1 = tan(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(nothomogeneous_exp_x)
{
	Var x = MakeVariable("x");

	auto f1 = exp(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}


BOOST_AUTO_TEST_CASE(nothomogeneous_sqrt_x)
{
	Var x = MakeVariable("x");

	auto f1 = sqrt(x);
	
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsPolynomial());
}



BOOST_AUTO_TEST_CASE(is_homogeneous_sin_0)
{
	Float n = MakeFloat("1");

	auto f1 = sin(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_cos_1)
{
	Float n = MakeFloat("1");

	auto f1 = cos(n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}

BOOST_AUTO_TEST_CASE(is_homogeneous_sin_1_plus_1)
{
	Float n = MakeFloat("1");

	auto f1 = sin(n + n);
	
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsPolynomial());
}




BOOST_AUTO_TEST_CASE(is_homogeneous_summands_homogeneous)
{
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * mpfr_float("4.12331") * pow(x,1);
	
	auto f1 = a+b;
	BOOST_CHECK(f1->IsHomogeneous());

	BOOST_CHECK(f1->IsHomogeneous(x));
	BOOST_CHECK(f1->IsHomogeneous(y));


	BOOST_CHECK(f1->IsPolynomial());

}

BOOST_AUTO_TEST_CASE(not_homogeneous_summands_inhomogeneous)
{
	Var x = MakeVariable("x");
	Var y = MakeVariable("y");

	auto a = pow(x,3) / 2;
	auto b = pow(x,2) * mpfr_float("4.12331");
	
	auto f1 = a+b;
	BOOST_CHECK(!f1->IsHomogeneous());

	BOOST_CHECK(!f1->IsHomogeneous(x));
	BOOST_CHECK( f1->IsHomogeneous(y));

	BOOST_CHECK(f1->IsPolynomial());
}


BOOST_AUTO_TEST_SUITE_END()





