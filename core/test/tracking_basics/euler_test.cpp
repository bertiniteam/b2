//This file is part of Bertini 2.
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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire





#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"

#include "trackers/ode_predictors.hpp"




extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(euler_predict_tracking_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;
using Float = bertini::node::Float;
using ExplicitRKPredictor = bertini::tracking::predict::ExplicitRKPredictor;

using Var = std::shared_ptr<Variable>;
using bertini::MakeVariable;
using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
using bertini::DefaultPrecision;


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
	Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	double norm_J, norm_J_inverse, size_proportion;
	
	Vec<dbl> predicted(2);
	predicted << dbl(2.40310963516214640018253210912048,0.187706567388887830930493342816564),
	dbl(0.370984337833979085688209698697074, 1.30889906180158745272421523674049);
	
	Vec<dbl> euler_prediction_result;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
	
	auto success_code = predictor->Predict(euler_prediction_result,
										   size_proportion,
										   norm_J, norm_J_inverse,
										   sys,
										   current_space, current_time,
										   delta_t,
										   condition_number_estimate,
										   num_steps_since_last_condition_number_computation,
										   frequency_of_CN_estimation,
										   tracking_tolerance,
										   AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
	for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
	BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	
	}
	
	
	BOOST_AUTO_TEST_CASE(circle_line_euler_mp)
	{
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
		
		// Starting time
		mpfr current_time("0.9");
		// Time step
		mpfr delta_t("-0.1");
		
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,2);
		AMP.coefficient_bound = 5;
		
		double norm_J, norm_J_inverse, size_proportion;
		
		Vec<mpfr> predicted(2);
		predicted << mpfr("2.40310963516214640018253210912048","0.187706567388887830930493342816564"),
		mpfr("0.370984337833979085688209698697074", "1.30889906180158745272421523674049");
		
		Vec<mpfr> euler_prediction_result;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(euler_prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_condition_number_computation,
											   frequency_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
		for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
			BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
		
	}
	
	
	
	
	
	
	
	BOOST_AUTO_TEST_CASE(monodromy_euler_d)
	{
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<dbl> current_space(2);
		current_space << dbl(4.641588833612776e-1), dbl(7.416198487095662e-1);
		
		// Starting time
		dbl current_time(0.7);
		// Time step
		dbl delta_t(-0.01);
		
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
		sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + mpfr_float("0.5")) );
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 2;
		
		
		Vec<dbl> predicted(2);
		predicted << dbl(0.417742995025149735840732480384),
		dbl(0.731506850772617663577442383933525);
		
		Vec<dbl> euler_prediction_result;
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(euler_prediction_result,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_condition_number_computation,
											   frequency_of_CN_estimation,
											   tracking_tolerance);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
		for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
			BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
		
		
	}
	
	
	
	
	BOOST_AUTO_TEST_CASE(monodromy_euler_mp)
	{
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
		
		// Starting time
		mpfr current_time("0.7");
		// Time step
		mpfr delta_t("-0.01");
		
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
		sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + mpfr_float("0.5")) );
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 2;
		
		
		Vec<mpfr> predicted(2);
		predicted << mpfr("0.417742995025149735840732480384"),
		mpfr("0.731506850772617663577442383934");
		
		Vec<mpfr> euler_prediction_result;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(euler_prediction_result,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_condition_number_computation,
											   frequency_of_CN_estimation,
											   tracking_tolerance);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
		for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
		 BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
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
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		AMP.coefficient_bound = 5;
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<dbl> prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::MatrixSolveFailureFirstPartOfPrediction);
		
		
	}
	
	
	
	BOOST_AUTO_TEST_CASE(euler_predict_linear_algebra_fails_mp)
	{
		// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
		
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("1.0"), mpfr("-4.0");
		
		// Starting time
		mpfr current_time(".7500000000000000000000000001");
		// Time step
		mpfr delta_t("-0.1");
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		AMP.coefficient_bound = 5;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<mpfr> prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::MatrixSolveFailureFirstPartOfPrediction);
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
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		AMP.coefficient_bound = 5;
		AMP.safety_digits_1 = 100;
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<dbl> prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_a_is_false_mp)
	{
		// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("1.0"), mpfr("-4.0");
		
		// Starting time
		mpfr current_time(".8");
		// Time step
		mpfr delta_t("-0.1");
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		AMP.coefficient_bound = 5;
		AMP.safety_digits_1 = 100;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<mpfr> prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
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
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		AMP.coefficient_bound = 5;
		AMP.safety_digits_2 = 100;
		
		AMP.SetPhiPsiFromBounds();
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<dbl> prediction_result;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	BOOST_AUTO_TEST_CASE(euler_predict_linear_criterion_c_is_false_mp)
	{
		// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
		
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("1.0"), mpfr("-4.0");
		
		// Starting time
		mpfr current_time(".8");
		// Time step
		mpfr delta_t("-0.1");
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x - 5*y) );
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		AMP.coefficient_bound = 5;
		AMP.safety_digits_2 = 100;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<mpfr> prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_cond_num_est,
											   freq_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
	
	
	BOOST_AUTO_TEST_CASE(circle_line_euler_change_precision)
	{
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
		
		// Starting time
		mpfr current_time("0.9");
		// Time step
		mpfr delta_t("-0.1");
		
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
		sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
		
		std::cout << "setting AMP config\n";

		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,2);
		AMP.coefficient_bound = 5;
		
		double norm_J, norm_J_inverse, size_proportion;
		
		Vec<mpfr> predicted(2);
		predicted << mpfr("2.40310963516214640018253210912048","0.187706567388887830930493342816564"),
		mpfr("0.370984337833979085688209698697074", "1.30889906180158745272421523674049");
		
		Vec<mpfr> euler_prediction_result;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(sys);
		
		auto success_code = predictor->Predict(euler_prediction_result,
											   size_proportion,
											   norm_J, norm_J_inverse,
											   sys,
											   current_space, current_time,
											   delta_t,
											   condition_number_estimate,
											   num_steps_since_last_condition_number_computation,
											   frequency_of_CN_estimation,
											   tracking_tolerance,
											   AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
		for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
			BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
		
		
		
		
		bertini::DefaultPrecision(50);
		Precision(current_space,50); assert(current_space(0).precision()==50);
		current_time.precision(50);
		delta_t.precision(50);
		sys.precision(50);assert(sys.precision()==50);
		Precision(predicted,50); assert(predicted(0).precision()==50);
		Precision(euler_prediction_result,50);  assert(euler_prediction_result(0).precision()==50);
		predictor->ChangePrecision(50);

		// Starting point in spacetime step
		current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
		
		// Starting time
		current_time = mpfr("0.9");
		// Time step
		delta_t = mpfr("-0.1");
		
		
		
		
		
		
		AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,2);
		AMP.coefficient_bound = 5;
		
		
		predicted << mpfr("2.4031096351621464001825321091204810500227008230702","0.18770656738888783093049334281656388282191769240875"),
		mpfr("0.37098433783397908568820969869707413571104273956141", "1.3088990618015874527242152367404908738825831867988");
		

		success_code = predictor->Predict(euler_prediction_result,
										  size_proportion,
										  norm_J, norm_J_inverse,
										  sys,
										  current_space, current_time,
										  delta_t,
										  condition_number_estimate,
										  num_steps_since_last_condition_number_computation,
										  frequency_of_CN_estimation,
										  tracking_tolerance,
										  AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(euler_prediction_result.size(),2);
		for (unsigned ii = 0; ii < euler_prediction_result.size(); ++ii)
		{
			BOOST_CHECK(abs(euler_prediction_result(ii)-predicted(ii)) < 1e-47);
		}
	}

BOOST_AUTO_TEST_SUITE_END()


























