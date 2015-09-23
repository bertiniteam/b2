//This file is part of Bertini 2.0.
//
//euler_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//euler_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with euler_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  euler_test.cpp
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

#include "tracking/predict.hpp"


using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = boost::multiprecision::mpfr_float;


template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;

extern double threshold_clearance_d;
extern boost::multiprecision::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(euler_predict_tracking_basics)


BOOST_AUTO_TEST_CASE(circle_line_euler_double)
{

	// Starting point in spacetime step    
	Vec<dbl> current_space(2);
	current_space << dbl(2.3,0.2), dbl(1.1, 1.87);

	// Starting time
	dbl current_time(0.9);
	// Time step
	dbl delta_t(-0.1);
	
	
	
	
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


	Vec<dbl> predicted(2);
	predicted << dbl(2.40310963516214640018253210912048,0.187706567388887830930493342816564),
		dbl(0.370984337833979085688209698697074, 1.30889906180158745272421523674049);

	Vec<dbl> euler_prediction_result;
	dbl next_time;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;

	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
								euler_prediction_result,
								sys,
								current_space, current_time, 
								delta_t,
								condition_number_estimate,
								num_steps_since_last_condition_number_computation, 
								frequency_of_CN_estimation, bertini::tracking::config::PrecisionType::Adaptive, 
								tracking_tolerance,
								AMP);
	
	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
	for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);

}


BOOST_AUTO_TEST_CASE(circle_line_euler_mp)
{
	boost::multiprecision::mpfr_float::default_precision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);

	// Starting point in spacetime step    
	Vec<mpfr> current_space(2);
	current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");

	// Starting time
	mpfr current_time("0.9");
	// Time step
	mpfr delta_t("-0.1");
	
	
	
	
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


	Vec<mpfr> predicted(2);
	predicted << mpfr("2.40310963516214640018253210912048","0.187706567388887830930493342816564"),
		mpfr("0.370984337833979085688209698697074", "1.30889906180158745272421523674049");

	Vec<mpfr> euler_prediction_result;
	mpfr next_time;
	
	mpfr_float tracking_tolerance("1e-5");
	mpfr_float condition_number_estimate;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;

	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
								euler_prediction_result,
								sys,
								current_space, current_time, 
								delta_t,
								condition_number_estimate,
								num_steps_since_last_condition_number_computation, 
								frequency_of_CN_estimation, bertini::tracking::config::PrecisionType::Adaptive, 
								tracking_tolerance,
								AMP);
	
	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
	for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);

}

	
	
	
	
BOOST_AUTO_TEST_CASE(monodromy_euler_mp)
{
	boost::multiprecision::mpfr_float::default_precision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("4.641588833612776e-1"), mpfr("7.416198487095662e-1");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1.0) + (1-t)*(pow(x,2) + 2.0) );
	sys.AddFunction( t*(pow(y,2)-1.0) + (1-t)*(pow(y,2) + .5) );
	
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 5;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.417742995025149735840732480383684"),
	mpfr("0.731506850772617663577442383933525");
	
	Vec<mpfr> euler_prediction_result;
	mpfr next_time;
	
	mpfr_float tracking_tolerance("1e-5");
	mpfr_float condition_number_estimate;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   euler_prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_condition_number_computation,
												   frequency_of_CN_estimation, bertini::tracking::config::PrecisionType::Adaptive,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code==bertini::tracking::SuccessCode::Success);
	BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
	for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	std::cout <<abs(euler_prediction_result(0)) << std::endl;
	
}

	
	
	



BOOST_AUTO_TEST_CASE(euler_predict_linear_algebra_fails_d)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	// Starting point in spacetime step
	Vec<dbl> current_space(2);
	current_space << dbl(1.0), dbl(-4.0);
	
	// Starting time
	dbl current_time(.75+1e-13);
	// Time step
	dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<dbl> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::Double,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::MatrixSolveFailure);
	

}
	
	

BOOST_AUTO_TEST_CASE(euler_predict_linear_algebra_fails_mp)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("1.0"), mpfr("-4.0");
	
	// Starting time
	mpfr current_time(".7500000000000000000000000001");
	// Time step
	mpfr delta_t("-0.1");
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	
	mpfr_float tracking_tolerance("1e-5");
	mpfr_float condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<mpfr> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::FixedMultiple,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::MatrixSolveFailure);
}


BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_a_is_false_d)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	// Starting point in spacetime step
	Vec<dbl> current_space(2);
	current_space << dbl(1.0), dbl(-4.0);
	
	// Starting time
	dbl current_time(.8);
	// Time step
	dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	AMP.safety_digits_1 = 100;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<dbl> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::Adaptive,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::HigherPrecisionNecessary);
}

BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_a_is_false_mp)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("1.0"), mpfr("-4.0");
	
	// Starting time
	mpfr current_time(".8");
	// Time step
	mpfr delta_t("-0.1");
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	AMP.safety_digits_1 = 100;
	
	mpfr_float tracking_tolerance("1e-5");
	mpfr_float condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<mpfr> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::Adaptive,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::HigherPrecisionNecessary);
}

BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_c_is_false_d)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	// Starting point in spacetime step
	Vec<dbl> current_space(2);
	current_space << dbl(1.0), dbl(-4.0);
	
	// Starting time
	dbl current_time(.8);
	// Time step
	dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	AMP.safety_digits_2 = 100;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<dbl> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::Adaptive,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::HigherPrecisionNecessary);
}

BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_c_is_false_mp)
{
	// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
	
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("1.0"), mpfr("-4.0");
	
	// Starting time
	mpfr current_time(".8");
	// Time step
	mpfr delta_t("-0.1");
	
	
	
	bertini::System sys;
	Var x = std::make_shared<Variable>("x"), y = std::make_shared<Variable>("y"), t = std::make_shared<Variable>("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
	
	auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
	
	AMP.coefficient_bound = 5;
	AMP.safety_digits_2 = 100;
	
	mpfr_float tracking_tolerance("1e-5");
	mpfr_float condition_number_estimate;
	
	unsigned num_steps_since_last_cond_num_est = 1;
	unsigned freq_of_CN_estimation = 1;
	
	Vec<mpfr> prediction_result;
	
	
	auto success_code = bertini::tracking::Predict(bertini::tracking::config::Predictor::Euler,
												   prediction_result,
												   sys,
												   current_space, current_time,
												   delta_t,
												   condition_number_estimate,
												   num_steps_since_last_cond_num_est,
												   freq_of_CN_estimation,
												   bertini::tracking::config::PrecisionType::Adaptive,
												   tracking_tolerance,
												   AMP);
	
	BOOST_CHECK(success_code == bertini::tracking::SuccessCode::HigherPrecisionNecessary);
}
BOOST_AUTO_TEST_SUITE_END()

























