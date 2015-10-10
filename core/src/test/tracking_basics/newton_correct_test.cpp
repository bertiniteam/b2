//This file is part of Bertini 2.0.
//
//newton_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//newton_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with newton_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015

//newton_test.cpp
//


#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"
#include "tracking/correct.hpp"


using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;

template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;


BOOST_AUTO_TEST_SUITE(newton_correct_tracking_basics)

BOOST_AUTO_TEST_CASE(circle_line_one_corrector_step_double)
{

	// Starting point in spacetime step    
	Vec<dbl> current_space(2);
	current_space << dbl(2.3,0.2), dbl(1.1, 1.87);

	// Starting time
	dbl current_time(0.9);
	// Time step
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;


	
	double tracking_tolerance(1e-5);


	Vec<dbl> corrected(2);
	corrected << dbl(1.36296628875178620892887063382866,0.135404746200380445814213878747082),
		dbl(0.448147673035459113010161338024478, -0.0193435351714829208306019826781546);

	Vec<dbl> newton_correction_result;

	tracking_tolerance = double(1e1);
	double path_truncation_threshold(1e4);
	unsigned max_num_newton_iterations = 1;
	unsigned min_num_newton_iterations = 1;
	auto success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               current_space, 
								               current_time, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               min_num_newton_iterations,
								               max_num_newton_iterations,
								               AMP);

	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
	for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_d);

}



BOOST_AUTO_TEST_CASE(circle_line_one_corrector_step_mp)
{

	// Starting point in spacetime step    
	Vec<mpfr> current_space(2);
	current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");

	// Starting time
	mpfr current_time("0.9");
	// Time step
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;


	
	bertini::mpfr_float tracking_tolerance("1e-5");


	Vec<mpfr> corrected(2);
	corrected << mpfr("1.36296628875178620892887063382866","0.135404746200380445814213878747082"),
		mpfr("0.448147673035459113010161338024478", "-0.0193435351714829208306019826781546");

	Vec<mpfr> newton_correction_result;

	tracking_tolerance = bertini::mpfr_float("1e1");
	bertini::mpfr_float path_truncation_threshold("1e4");
	unsigned max_num_newton_iterations = 1;
	unsigned min_num_newton_iterations = 1;
	auto success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               current_space, 
								               current_time, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               min_num_newton_iterations,
								               max_num_newton_iterations,
								               AMP);

	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
	for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_mp);

}



BOOST_AUTO_TEST_CASE(circle_line_two_corrector_steps_double)
{

	// Starting point in spacetime step    
	Vec<dbl> current_space(2);
	current_space << dbl(2.3,0.2), dbl(1.1, 1.87);

	// Starting time
	dbl current_time(0.9);
	// Time step
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;


	
	double tracking_tolerance(1e-5);


	Vec<dbl> corrected(2);
	corrected << dbl(1.14542104415948767661671388923986, 0.0217584797792294631577151109622764),
		dbl(0.47922556512007318905475515868002, -0.00310835425417563759395930156603948);

	Vec<dbl> newton_correction_result;

	tracking_tolerance = double(1e1);
	double path_truncation_threshold(1e4);
	unsigned max_num_newton_iterations = 2;
	unsigned min_num_newton_iterations = 2;
	auto success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               current_space, 
								               current_time, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               min_num_newton_iterations,
								               max_num_newton_iterations,
								               AMP);

	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
	for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_d);

}



BOOST_AUTO_TEST_CASE(circle_line_two_corrector_steps_mp)
{

	// Starting point in spacetime step    
	Vec<mpfr> current_space(2);
	current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");

	// Starting time
	mpfr current_time("0.9");
	// Time step
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	

	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);

	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;


	
	bertini::mpfr_float tracking_tolerance("1e-5");


	Vec<mpfr> corrected(2);
	corrected << mpfr("1.14542104415948767661671388923986", "0.0217584797792294631577151109622764"),
		mpfr("0.47922556512007318905475515868002", "-0.00310835425417563759395930156603948");

	Vec<mpfr> newton_correction_result;

	tracking_tolerance = bertini::mpfr_float("1e1");
	bertini::mpfr_float path_truncation_threshold("1e4");
	unsigned max_num_newton_iterations = 2;
	unsigned min_num_newton_iterations = 2;
	auto success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               current_space, 
								               current_time, 
								               tracking_tolerance,
								               path_truncation_threshold,
								               min_num_newton_iterations,
								               max_num_newton_iterations,
								               AMP);

	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
	for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_mp);

}






BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_B_violated_double)
{
	BOOST_CHECK("implemented case where newton step requests higher precision due to AMP criterion B"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_B_violated_mp)
{
	BOOST_CHECK("implemented case where newton step requests higher precision due to AMP criterion B"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_C_violated_double)
{
	BOOST_CHECK("implemented case where newton step requests higher precision due to AMP criterion C"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_C_violated_mp)
{
	BOOST_CHECK("implemented case where newton step requests higher precision due to AMP criterion C"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_linear_algebra_fails_double)
{
	BOOST_CHECK("implemented case where newton step linear algebra fails"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_linear_algebra_fails_mp)
{
	BOOST_CHECK("implemented case where newton step linear algebra fails"=="true");
}


BOOST_AUTO_TEST_CASE(newton_step_going_to_infinity_d)
{
	BOOST_CHECK("implemented case where newton loop terminates due to going to infinity"=="true");
}

BOOST_AUTO_TEST_CASE(newton_step_going_to_infinity_mp)
{
	BOOST_CHECK("implemented case where newton loop terminates due to going to infinity"=="true");
}


BOOST_AUTO_TEST_SUITE_END()





