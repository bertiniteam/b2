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

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"
#include "tracking.hpp"


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




BOOST_AUTO_TEST_SUITE(tracking_basics)

BOOST_AUTO_TEST_CASE(circle_line_euler_one_corr_step_double)
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
	sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
	

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
	unsigned num_steps_since_last_condition_number_computation = 2;
	unsigned frequency_of_CN_estimation = 2;

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


	next_time = current_time + delta_t;

	//////////////
	//
	//  begin correction
	//
	////////


	Vec<dbl> corrected(2);
	corrected << dbl(1.46335169441695765372296267763834,0.101522740428208548465501963605153),
		dbl(0.119255179018453854728230516080369, -0.0225606089840463441034448808011452);

	Vec<dbl> newton_correction_result;

	tracking_tolerance = double(1e1);
	double path_truncation_threshold(1e4);
	unsigned max_num_newton_iterations = 1;
	unsigned min_num_newton_iterations = 1;
	success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               euler_prediction_result, 
								               next_time, 
								               bertini::tracking::config::PrecisionType::Adaptive, 
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





BOOST_AUTO_TEST_CASE(circle_line_euler_one_corr_step_mp)
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
	sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
	

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


	next_time = current_time + delta_t;

	//////////////
	//
	//  begin correction
	//
	////////


	Vec<mpfr> corrected(2);
	corrected << mpfr("1.46335169441695765372296267763834","0.101522740428208548465501963605153"),
		mpfr("0.119255179018453854728230516080369", "-0.0225606089840463441034448808011452");

	Vec<mpfr> newton_correction_result;

	tracking_tolerance = mpfr_float("1e1");
	mpfr_float path_truncation_threshold("1e4");
	unsigned max_num_newton_iterations = 1;
	unsigned min_num_newton_iterations = 1;
	success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               euler_prediction_result, 
								               next_time, 
								               bertini::tracking::config::PrecisionType::Adaptive, 
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









BOOST_AUTO_TEST_CASE(circle_line_euler_two_corr_steps_double)
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
	sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
	

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


	next_time = current_time + delta_t;

	//////////////
	//
	//  begin correction
	//
	////////


	Vec<dbl> corrected(2);
	corrected << dbl(1.27408462511795319743443276342357,0.0133086443984782567063605708698885),
		dbl(0.161314527751565956125681608128097, -0.0029574765329951681569690157488641);

	Vec<dbl> newton_correction_result;

	tracking_tolerance = double(1e1);
	double path_truncation_threshold(1e4);
	unsigned max_num_newton_iterations = 2;
	unsigned min_num_newton_iterations = 2;
	success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               euler_prediction_result, 
								               next_time, 
								               bertini::tracking::config::PrecisionType::Adaptive, 
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


BOOST_AUTO_TEST_CASE(circle_line_euler_two_corr_steps_mp)
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
	sys.AddFunction( t*(y-1.0) + (1-t)*(2.0*x + 5.0*y) );
	

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


	next_time = current_time + delta_t;

	//////////////
	//
	//  begin correction
	//
	////////


	Vec<mpfr> corrected(2);
	corrected << mpfr("1.27408462511795319743443276342357","0.0133086443984782567063605708698885"),
		mpfr("0.161314527751565956125681608128097", "-0.0029574765329951681569690157488641");

	Vec<mpfr> newton_correction_result;

	tracking_tolerance = mpfr_float("1e1");
	mpfr_float path_truncation_threshold("1e4");
	unsigned max_num_newton_iterations = 2;
	unsigned min_num_newton_iterations = 2;
	success_code = bertini::tracking::Correct(newton_correction_result,
								               sys,
								               euler_prediction_result, 
								               next_time, 
								               bertini::tracking::config::PrecisionType::Adaptive, 
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






BOOST_AUTO_TEST_SUITE_END()

























