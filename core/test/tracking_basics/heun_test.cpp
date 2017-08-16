//This file is part of Bertini 2.0.
//
//heun_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//heun_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with heun_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  Heun_test.cpp
//
//  copyright 2016
//  James B. Collins
//  West Texas A&M University
//  Department of Mathematics
//  Spring 2016



#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "limbo.hpp"
#include "mpfr_complex.hpp"

#include "trackers/ode_predictors.hpp"




extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(heun_predict_tracking_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;
using Float = bertini::node::Float;
using ExplicitRKPredictor = bertini::tracking::predict::ExplicitRKPredictor;

using bertini::MakeVariable;
using bertini::MakeFloat;
using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;


BOOST_AUTO_TEST_CASE(circle_line_heun_double)
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
	
	double norm_J, norm_J_inverse, size_proportion, error_est;
	
	Vec<dbl> predicted(2);
	predicted << dbl(2.38948874619536140814029774733947,0.208678935223681033727262214382917),
	dbl(0.524558056401030798191044945035673, 1.43029356995029310361616395235936);
	double predicted_error = .197349645229023708608160063982175;
	
	Vec<dbl> heun_prediction_result;
	dbl next_time;
	
	double tracking_tolerance(1e-5);
	double condition_number_estimate;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
	
	auto success_code = predictor->Predict(heun_prediction_result,
										   error_est,
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
	BOOST_CHECK_EQUAL(heun_prediction_result.size(),2);
	for (unsigned ii = 0; ii < heun_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(heun_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	BOOST_CHECK(fabs(error_est / predicted_error - 1) < threshold_clearance_d);
	
	}
	
	
	
	
	
	
	
	
	
	
	BOOST_AUTO_TEST_CASE(circle_line_heun_mp)
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
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		Vec<mpfr> predicted(2);
		predicted << mpfr("2.38948874619536140814029774733947","0.208678935223681033727262214382917"),
		mpfr("0.524558056401030798191044945035673", "1.43029356995029310361616395235936");
		double predicted_error = double(.197349645229023708608160063982175);
		
		Vec<mpfr> heun_prediction_result;
		mpfr next_time;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		BOOST_CHECK_EQUAL(heun_prediction_result.size(),2);
		for (unsigned ii = 0; ii < heun_prediction_result.size(); ++ii)
			BOOST_CHECK(abs(heun_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
		
		using std::abs;
		BOOST_CHECK(abs(error_est - predicted_error) < std::numeric_limits<double>::epsilon());	
	}
	
	
	
	
	
	
	
	BOOST_AUTO_TEST_CASE(monodromy_heun_d)
	{
		bertini::DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		
		// Starting point in spacetime step
		Vec<dbl> current_space(2);
		current_space << dbl(0.464158883361277585510862309093), dbl(0.74161984870956629487113974408);
		
		// Starting time
		dbl current_time(0.7);
		// Time step
		dbl delta_t(-0.01);
		
		
		
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		std::shared_ptr<Float> half = MakeFloat("0.5");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
		sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 2;
		
		
		Vec<dbl> predicted(2);
		predicted << dbl(0.412299156269677938503694812160886),
		dbl(0.731436945256924470273568899877140);
		double predicted_error = 0.00544428757292458409463632380167773;
		
		Vec<dbl> heun_prediction_result;
		double next_time;
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		BOOST_CHECK_EQUAL(heun_prediction_result.size(),2);
		for (unsigned ii = 0; ii < heun_prediction_result.size(); ++ii)
		{
			BOOST_CHECK(abs(heun_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
		}
		
		BOOST_CHECK(fabs(error_est / predicted_error - 1) < threshold_clearance_d);
		
	}
	
	
	
	BOOST_AUTO_TEST_CASE(monodromy_heun_mp)
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
		std::shared_ptr<Float> half = MakeFloat("0.5");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		// Define homotopy system
		sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
		sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 2;
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		
		Vec<mpfr> predicted(2);
		predicted << mpfr("0.412299156269677938503694812160886"),
		mpfr("0.731436945256924470273568899877140");
		double predicted_error = double(0.00544428757292458409463632380167773);
		
		Vec<mpfr> heun_prediction_result;
		mpfr next_time;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		BOOST_CHECK_EQUAL(heun_prediction_result.size(),2);
		for (unsigned ii = 0; ii < heun_prediction_result.size(); ++ii)
		{
			BOOST_CHECK(abs(heun_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
		}
		
		using std::abs;
		BOOST_CHECK(abs(error_est / predicted_error - 1) < std::numeric_limits<double>::epsilon());
	}
	
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_algebra_fails_d)
	{
		// Circle line homotopy has singular point at (x,y) = (1,-4) and t = .75
		
		// Starting point in spacetime step
		Vec<dbl> current_space(2);
		current_space << dbl(1.0), dbl(-4.0);
		
		// Starting time
		dbl current_time(.75);
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
		
		double tracking_tolerance(1e-5);
		double condition_number_estimate;
		
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		Vec<dbl> heun_prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::MatrixSolveFailureFirstPartOfPrediction);
		
	}
	
	
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_algebra_fails_mp)
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
		
		double norm_J, norm_J_inverse, size_proportion, error_est;
		
		AMP.coefficient_bound = 5;
		
		double tracking_tolerance = 1e-5;
		double condition_number_estimate;
		
		unsigned num_steps_since_last_cond_num_est = 1;
		unsigned freq_of_CN_estimation = 1;
		
		Vec<mpfr> heun_prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::MatrixSolveFailureFirstPartOfPrediction);
	}
	
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_criterion_a_is_false_d)
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
		
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		Vec<dbl> heun_prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_criterion_a_is_false_mp)
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
		
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		Vec<mpfr> heun_prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_criterion_c_is_false_d)
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
		
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		Vec<dbl> heun_prediction_result;
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	BOOST_AUTO_TEST_CASE(heun_predict_linear_criterion_c_is_false_mp)
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
		
		unsigned num_steps_since_last_condition_number_computation = 1;
		unsigned frequency_of_CN_estimation = 1;
		
		Vec<mpfr> heun_prediction_result;
		
		
		std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::HeunEuler,sys);
		
		auto success_code = predictor->Predict(heun_prediction_result,
											   error_est,
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
		
		BOOST_CHECK(success_code == bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
BOOST_AUTO_TEST_SUITE_END()


























