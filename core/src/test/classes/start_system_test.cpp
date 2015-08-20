//This file is part of Bertini 2.0.
//
//start_system_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  start_system_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

//start_system_test.cpp
//


#include <boost/test/unit_test.hpp>



#include "start_system.hpp"

using System = bertini::System;
using Var = std::shared_ptr<bertini::Variable>;
using VariableGroup = bertini::VariableGroup;


template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

double threshold_clearance_d = 1e-15;
double threshold_clearance_mp = 1e-27;
unsigned FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS = 30;

BOOST_AUTO_TEST_SUITE(system_class)





BOOST_AUTO_TEST_CASE(make_total_degree_system_linear)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");

	VariableGroup v;
	v.push_back(x); v.push_back(y);

	sys.AddVariableGroup(v);
	sys.AddFunction(x + y - 1);
	sys.AddFunction(x - 0.5*y - 1);


	bertini::start_system::TotalDegree TD(sys);

	auto d = TD.Degrees();

	BOOST_CHECK_EQUAL(d.size(),2);
	if (d.size()==2)
	{
		BOOST_CHECK_EQUAL(d[0],1);
		BOOST_CHECK_EQUAL(d[1],1);
	}
	
	
	BOOST_CHECK_EQUAL(TD.NumVariables(),2);

	BOOST_CHECK_EQUAL(TD.NumStartPoints(), 1);
}



BOOST_AUTO_TEST_CASE(make_total_degree_system_quadratic)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");

	VariableGroup v;
	v.push_back(x); v.push_back(y);

	sys.AddVariableGroup(v);
	sys.AddFunction(x*y + y - 1);
	sys.AddFunction(x*x - 0.5*y - x*y);


	bertini::start_system::TotalDegree TD(sys);

	auto d = TD.Degrees();

	BOOST_CHECK_EQUAL(d.size(),2);
	if (d.size()==2)
	{
		BOOST_CHECK_EQUAL(d[0],2);
		BOOST_CHECK_EQUAL(d[1],2);
	}
	
	
	BOOST_CHECK_EQUAL(TD.NumVariables(),2);
	BOOST_CHECK_EQUAL(TD.NumStartPoints(), 4);
}



BOOST_AUTO_TEST_CASE(linear_total_degree_start_system)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+1);
	sys.AddFunction(x+y+bertini::Pi());

	bertini::start_system::TotalDegree TD(sys);

	auto deg = TD.Degrees();

	BOOST_CHECK_EQUAL(deg.size(),2);
	if (deg.size()==2)
	{
		BOOST_CHECK_EQUAL(deg[0],1);
		BOOST_CHECK_EQUAL(deg[1],1);
	}

	Vec<dbl> vals(2);
	vals << dbl(1.0),dbl(1.0);

	auto sysvals = TD.Eval(vals);

	for (unsigned ii = 0; ii < 2; ++ii)
		BOOST_CHECK( abs(sysvals(ii) - (1.0 - dbl(TD.RandomValue(ii)))) < threshold_clearance_d);



	auto J = TD.Jacobian(vals);

	BOOST_CHECK_EQUAL(J(0,0),1.0);
	BOOST_CHECK_EQUAL(J(0,1),0.0);

	BOOST_CHECK_EQUAL(J(1,0),0.0);
	BOOST_CHECK_EQUAL(J(1,1),1.0);


	vals << dbl(0.0),dbl(0.0);

	sysvals = TD.Eval(vals);

	for (unsigned ii = 0; ii < 2; ++ii)
		BOOST_CHECK(abs(sysvals(ii)+dbl(TD.RandomValue(ii))) < threshold_clearance_d);
	

	J = TD.Jacobian(vals);

	BOOST_CHECK_EQUAL(J(0,0),1.0);
	BOOST_CHECK_EQUAL(J(0,1),0.0);

	BOOST_CHECK_EQUAL(J(1,0),0.0);
	BOOST_CHECK_EQUAL(J(1,1),1.0);

	BOOST_CHECK_EQUAL(TD.NumStartPoints(), 1);
}







BOOST_AUTO_TEST_CASE(quadratic_cubic_quartic_total_degree_start_system)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y"), z = std::make_shared<bertini::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + 0.5);
	sys.AddFunction(pow(x,3)+x*y+bertini::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	bertini::start_system::TotalDegree TD(sys);

	auto deg = TD.Degrees();

	BOOST_CHECK_EQUAL(deg.size(),3);
	if (deg.size()==3)
	{
		BOOST_CHECK_EQUAL(deg[0],2);
		BOOST_CHECK_EQUAL(deg[1],3);
		BOOST_CHECK_EQUAL(deg[2],4);
	}

	Vec<dbl> vals(3);
	vals << dbl(1.0),dbl(1.0),dbl(1.0);

	auto sysvals = TD.Eval(vals);

	for (unsigned ii = 0; ii < 3; ++ii)
		BOOST_CHECK( abs(sysvals(ii) - (1.0 - dbl(TD.RandomValue(ii)))) < threshold_clearance_d);


	auto J = TD.Jacobian(vals);

	BOOST_CHECK_EQUAL(J(0,0),2.0);
	BOOST_CHECK_EQUAL(J(0,1),0.0);
	BOOST_CHECK_EQUAL(J(0,2),0.0);

	BOOST_CHECK_EQUAL(J(1,0),0.0);
	BOOST_CHECK_EQUAL(J(1,1),3.0);
	BOOST_CHECK_EQUAL(J(1,2),0.0);

	BOOST_CHECK_EQUAL(J(2,0),0.0);
	BOOST_CHECK_EQUAL(J(2,1),0.0);
	BOOST_CHECK_EQUAL(J(2,2),4.0);



	vals << dbl(0.0),dbl(0.0),dbl(0.0);

	sysvals = TD.Eval(vals);

	for (unsigned ii = 0; ii < 3; ++ii)
		BOOST_CHECK(abs(sysvals(ii)+dbl(TD.RandomValue(ii))) < threshold_clearance_d);

	J = TD.Jacobian(vals);

	BOOST_CHECK_EQUAL(J(0,0),0.0);
	BOOST_CHECK_EQUAL(J(0,1),0.0);
	BOOST_CHECK_EQUAL(J(0,2),0.0);

	BOOST_CHECK_EQUAL(J(1,0),0.0);
	BOOST_CHECK_EQUAL(J(1,1),0.0);
	BOOST_CHECK_EQUAL(J(1,2),0.0);

	BOOST_CHECK_EQUAL(J(2,0),0.0);
	BOOST_CHECK_EQUAL(J(2,1),0.0);
	BOOST_CHECK_EQUAL(J(2,2),0.0);


	BOOST_CHECK_EQUAL(TD.NumStartPoints(), 24);
}




BOOST_AUTO_TEST_CASE(quadratic_cubic_quartic_start_points)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y"), z = std::make_shared<bertini::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + 0.5);
	sys.AddFunction(pow(x,3)+x*y+bertini::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	bertini::start_system::TotalDegree TD(sys);

	for (size_t ii = 0; ii < TD.NumStartPoints(); ++ii)
	{
		auto start = TD.StartPoint<dbl>(ii);
		auto function_values = TD.Eval(start);

		for (size_t jj = 0; jj < function_values.size(); ++jj)
			BOOST_CHECK(abs(function_values(jj)) < threshold_clearance_d);
	}

	boost::multiprecision::mpfr_float::default_precision(FUNCTION_TREE_TEST_MPFR_DEFAULT_DIGITS);

	for (size_t ii = 0; ii < TD.NumStartPoints(); ++ii)
	{
		auto start = TD.StartPoint<mpfr>(ii);
		auto function_values = TD.Eval(start);

		for (size_t jj = 0; jj < function_values.size(); ++jj)
		{
			BOOST_CHECK(abs(function_values(jj)) < threshold_clearance_mp);
		}
	}

}



BOOST_AUTO_TEST_CASE(quadratic_cubic_quartic_all_the_way_to_final_system)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::Variable>("x"), y = std::make_shared<bertini::Variable>("y"), z = std::make_shared<bertini::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + 0.5);
	sys.AddFunction(pow(x,3)+x*y+bertini::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	bertini::start_system::TotalDegree TD(sys);

	Var t = std::make_shared<bertini::Variable>("t");

	auto final_mixed_sum = (1-t) * sys + t * TD;
	final_mixed_sum.AddPathVariable(t);
	final_mixed_sum.Homogenize();

	BOOST_CHECK_EQUAL(final_mixed_sum.NumVariables(),4);
	BOOST_CHECK(final_mixed_sum.IsHomogeneous());
	BOOST_CHECK(final_mixed_sum.IsPolynomial());

}




BOOST_AUTO_TEST_SUITE_END()




