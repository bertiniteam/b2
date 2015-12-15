//This file is part of Bertini 2.0.
//
//amp_criteria_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_criteria_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_criteria_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  amp_criteria_test.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015



#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"

#include "tracking/amp_criteria.hpp"


using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using mpq_rational = bertini::mpq_rational;
using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(amp_criteria_tracking_basics)


BOOST_AUTO_TEST_CASE(AMP_criteriaA_double)
{
	/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<dbl> current_space(2);
	current_space << dbl(256185069753.4088,-387520022558.0519),
					 dbl(-0.021,-0.177);

	dbl current_time(0,0);
	dbl delta_t(.1,0);
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	Mat<dbl> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();

	Vec<dbl> randy = Vec<dbl>::Random(sys.NumVariables());
	Vec<dbl> temp_soln = LU.solve(randy);
					
	auto norm_J = dh_dx.norm();
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion A failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_1 = 32000;

	auto CritA = bertini::tracking::amp::CriterionA(norm_J,norm_J_inverse,AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritA,false);
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaA_mp)
{
		/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<mpfr> current_space(2);
	current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
					 mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");

	mpfr current_time("0");
	mpfr delta_t(".1");
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	Mat<mpfr> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();

	Vec<mpfr> randy = Vec<mpfr>::Random(sys.NumVariables());
	Vec<mpfr> temp_soln = LU.solve(randy);
					
	auto norm_J = dh_dx.norm();
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion A failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_1 = 32000;

	auto CritA = bertini::tracking::amp::CriterionA(norm_J,norm_J_inverse,AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritA,false);

}


BOOST_AUTO_TEST_CASE(AMP_criteriaB_double)
{
		/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<dbl> current_space(2);
	current_space << dbl(256185069753.4088,-387520022558.0519),
					 dbl(-0.021,-0.177);

	dbl current_time(0,0);
	dbl delta_t(.1,0);
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	auto f = sys.Eval(current_space, current_time);
	Mat<dbl> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();
	auto delta_z = LU.solve(-f);

	Vec<dbl> randy = Vec<dbl>::Random(sys.NumVariables());
	Vec<dbl> temp_soln = LU.solve(randy);
					
	auto norm_J = dh_dx.norm();
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion B failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_1 = 32000;

	unsigned int num_newton_iterations_remaining = 1;
	auto TrackTolBeforeEG = 10e-5; //Obtained from Bertini Book.

	auto CritB = bertini::tracking::amp::CriterionB(norm_J,norm_J_inverse,num_newton_iterations_remaining,TrackTolBeforeEG,delta_z.norm(),AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritB,false);
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaB_mp)
{
	/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<mpfr> current_space(2);
	current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
					 mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");

	mpfr current_time("0");
	mpfr delta_t(".1");
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	auto f = sys.Eval(current_space, current_time);
	Mat<mpfr> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();
	auto delta_z = LU.solve(-f);

	Vec<mpfr> randy = Vec<mpfr>::Random(sys.NumVariables());
	Vec<mpfr> temp_soln = LU.solve(randy);
					
	auto norm_J = dh_dx.norm();
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion B failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_1 = 32000;
	unsigned int num_newton_iterations_remaining = 1;
	mpfr_float TrackTolBeforeEG = mpfr_float("10e-5"); //Obtained from Bertini Book.

	auto CritB = bertini::tracking::amp::CriterionB(norm_J,norm_J_inverse,num_newton_iterations_remaining,TrackTolBeforeEG,delta_z.norm(),AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritB,false);
}

BOOST_AUTO_TEST_CASE(AMP_criteriaC_double)
{
	/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<dbl> current_space(2);
	current_space << dbl(256185069753.4088,-387520022558.0519),
					 dbl(-0.021,-0.177);

	dbl current_time(0,0);
	dbl delta_t(.1,0);
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	Mat<dbl> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();

	Vec<dbl> randy = Vec<dbl>::Random(sys.NumVariables());
	Vec<dbl> temp_soln = LU.solve(randy);
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion C failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_2 = 32000;
	auto TrackTolBeforeEG = 10e-5; //Obtained from Bertini Book.

	auto CritC = bertini::tracking::amp::CriterionC(norm_J_inverse,current_space,TrackTolBeforeEG,AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritC,false);
}
	
BOOST_AUTO_TEST_CASE(AMP_criteriaC_mp)
{
		/*
	Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict 
	to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. 
	Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition. 
	*/
	//Setting upt current space and time values for evaluation
	Vec<mpfr> current_space(2);
	current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
					 mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");

	mpfr current_time("0");
	mpfr delta_t(".1");
	current_time += delta_t;

	//Defining the system and variables. 
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	VariableGroup vars{x,y};

	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
	sys.AddFunction(y - pow(x,2));

	//For Criterion A to be checked we need Norm_J and inverse of Norm_J these were taken from Euler.hpp
	Mat<mpfr> dh_dx = sys.Jacobian(current_space, current_time); 
	auto LU = dh_dx.lu();

	Vec<mpfr> randy = Vec<mpfr>::Random(sys.NumVariables());
	Vec<mpfr> temp_soln = LU.solve(randy);
	auto norm_J_inverse = temp_soln.norm();


	//Setting up saftety digits to trigger AMP Criterion B failure.
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	AMP.safety_digits_2 = 32000;
	mpfr_float TrackTolBeforeEG = mpfr_float("10e-5"); //Obtained from Bertini Book.

	auto CritC = bertini::tracking::amp::CriterionC(norm_J_inverse,current_space,TrackTolBeforeEG,AMP);


	//Check to make sure we failed.
	BOOST_CHECK_EQUAL(CritC,false);
}





BOOST_AUTO_TEST_SUITE_END()

