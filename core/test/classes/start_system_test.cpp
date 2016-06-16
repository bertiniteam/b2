//This file is part of Bertini 2.
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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


#include <boost/test/unit_test.hpp>



#include "bertini2/start_system.hpp"

using System = bertini::System;

using Variable = bertini::node::Variable;
using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;

using mpq_rational = bertini::mpq_rational;
using mpfr_float = bertini::mpfr_float;
using mpz_int = bertini::mpz_int;
using dbl = bertini::dbl;
using mpfr = bertini::mpfr;

template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
extern double relaxed_threshold_clearance_d;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned CLASS_TEST_MPFR_DEFAULT_DIGITS;


BOOST_AUTO_TEST_SUITE(system_class)





BOOST_AUTO_TEST_CASE(make_total_degree_system_linear)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y");

	VariableGroup v;
	v.push_back(x); v.push_back(y);

	sys.AddVariableGroup(v);
	sys.AddFunction(x + y - 1);
	sys.AddFunction(x - mpfr_float("0.5")*y - 1);


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
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y");

	VariableGroup v;
	v.push_back(x); v.push_back(y);

	sys.AddVariableGroup(v);
	sys.AddFunction(x*y + y - 1);
	sys.AddFunction(x*x - mpfr_float("0.5")*y - x*y);


	bertini::start_system::TotalDegree TD(sys);

	auto d = TD.Degrees();

	BOOST_CHECK_EQUAL(TD.NumVariables(),2);
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
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y");

	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+1);
	sys.AddFunction(x+y+bertini::node::Pi());

	bertini::start_system::TotalDegree TD(sys);

	auto deg = TD.Degrees();

	BOOST_CHECK_EQUAL(TD.NumVariables(),2);
	BOOST_CHECK(!TD.IsPatched());

	BOOST_CHECK_EQUAL(deg.size(),2);
	
	BOOST_CHECK_EQUAL(deg[0],1);
	BOOST_CHECK_EQUAL(deg[1],1);
	
	VariableGroup variable_ordering = TD.Variables();

	BOOST_CHECK_EQUAL(variable_ordering.size(), 2);

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
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
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
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	bertini::start_system::TotalDegree TD(sys);

	for (mpz_int ii = 0; ii < TD.NumStartPoints(); ++ii)
	{
		auto start = TD.StartPoint<dbl>(ii);
		auto function_values = TD.Eval(start);
		
		for (size_t jj = 0; jj < function_values.size(); ++jj)
			BOOST_CHECK(abs(function_values(jj)) < relaxed_threshold_clearance_d);
	}

	bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	for (mpz_int ii = 0; ii < TD.NumStartPoints(); ++ii)
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
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars{x,y,z};

	sys.AddVariableGroup(vars);  
	sys.AddFunction(y+x*y + mpfr_float("0.5"));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	bertini::start_system::TotalDegree TD(sys);

	Var t = std::make_shared<bertini::node::Variable>("t");

	auto final_mixed_sum = (1-t) * sys + t * TD;
	final_mixed_sum.AddPathVariable(t);

	BOOST_CHECK(!final_mixed_sum.IsHomogeneous());
	BOOST_CHECK(!final_mixed_sum.IsPatched());

	final_mixed_sum.Homogenize();
	final_mixed_sum.AutoPatch();

	BOOST_CHECK_EQUAL(final_mixed_sum.NumVariables(),4);
	BOOST_CHECK_EQUAL(final_mixed_sum.NumNaturalVariables(),3);
	BOOST_CHECK_EQUAL(final_mixed_sum.NumFunctions(),3);
	BOOST_CHECK_EQUAL(final_mixed_sum.NumTotalFunctions(),4);

	BOOST_CHECK(final_mixed_sum.IsHomogeneous());
	BOOST_CHECK(final_mixed_sum.IsPolynomial());
	BOOST_CHECK(final_mixed_sum.IsPatched());

	Vec<mpfr> v(4);
	v << mpfr(1), mpfr(1), mpfr(1), mpfr(1);

	auto f = final_mixed_sum.Eval(v,mpfr::rand());
	BOOST_CHECK_EQUAL(f.size(), 4);

	auto J = final_mixed_sum.Jacobian(v,mpfr::rand());
	BOOST_CHECK_EQUAL(J.rows(), 4);
	BOOST_CHECK_EQUAL(J.cols(), 4);	
}





BOOST_AUTO_TEST_CASE(start_system_total_degree_nonpolynomial_should_throw)
{
	bertini::System sys;
	Var x = std::make_shared<bertini::node::Variable>("x"), y = std::make_shared<bertini::node::Variable>("y"), z = std::make_shared<bertini::node::Variable>("z");

	VariableGroup vars;
	vars.push_back(x); vars.push_back(y); vars.push_back(z);

	sys.AddVariableGroup(vars);  
	sys.AddFunction(exp(y)+x*y + mpq_rational(1,2));
	sys.AddFunction(pow(x,3)+x*y+bertini::node::E());
	sys.AddFunction(pow(x,2)*pow(y,2)+x*y*z*z - 1);

	BOOST_CHECK_THROW(bertini::start_system::TotalDegree TD(sys), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(total_degree_start_system_coefficient_bound_degree_bound)
{
	/*
	In this example we take a decoupled system, homogenize and patch it. Track to endgame boundary and then run our endgame on the space
	values we have. 

	Take note that in this example we have two successes while the other hit a min track time. This is accounted for. 
	*/
	DefaultPrecision(30);

	Var x = std::make_shared<Variable>("x");
	Var y = std::make_shared<Variable>("y");
	Var t = std::make_shared<Variable>("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddVariableGroup(v);

	sys.AddFunction(pow(x-1,3));
	sys.AddFunction(pow(y-1,2));
	sys.Homogenize();
	sys.AutoPatch();

	BOOST_CHECK(sys.IsHomogeneous());
	BOOST_CHECK(sys.IsPatched());	

	

	auto TD = bertini::start_system::TotalDegree(sys);
	TD.Homogenize();
	BOOST_CHECK(TD.IsHomogeneous());
	BOOST_CHECK(TD.IsPatched());


	auto final_system = (1-t)*sys + t*TD;
	final_system.AddPathVariable(t);
}



BOOST_AUTO_TEST_SUITE_END()




