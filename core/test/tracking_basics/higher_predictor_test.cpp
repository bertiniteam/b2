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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire





#include <boost/test/unit_test.hpp>

#include <boost/multiprecision/mpfr.hpp>
#include "bertini2/mpfr_complex.hpp"

#include "bertini2/trackers/ode_predictors.hpp"




extern double threshold_clearance_d;
extern bertini::real_mp threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(higher_predict_tracking_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;
using Complex = bertini::node::Complex;
using ExplicitRKPredictor = bertini::tracking::predict::ExplicitRKPredictor;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using complex_dbl = std::complex<double>;
using mpfr = bertini::complex_mp;
using real_mp = bertini::real_mp;


template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;


using bertini::Precision;
using bertini::DefaultPrecision;


template<typename NumT, typename ...T>
NumT NumFromString(T... s)
{return bertini::NumTraits<NumT>::FromString(s...);}

using std::abs;

using NumErrorT = bertini::NumErrorT;


//////////////////////////////////////////////
//
//	RK4
//
////////////////////////
BOOST_AUTO_TEST_CASE(circle_line_RK4_double)
{
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(2.3,0.2), complex_dbl(1.1, 1.87);
	
	// Starting time
	complex_dbl current_time(0.9);
	// Time step
	complex_dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(2.39187197874999601460772208561997,0.215631510575697758920211277830812),
	complex_dbl(0.524028449166667552309395092084449, 1.42874855320540049801773082714871);
	
	Vec<complex_dbl> RK4_prediction_result;

	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RK4,sys);
	
	auto success_code = predictor->Predict(RK4_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RK4_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RK4_prediction_result.size(); ++ii)
	BOOST_CHECK(abs(RK4_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	
	}


BOOST_AUTO_TEST_CASE(circle_line_RK4_mp)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39187197874999601460772208561997","0.215631510575697758920211277830812"),
	mpfr("0.524028449166667552309395092084449", "1.42874855320540049801773082714871");
	
	Vec<mpfr> RK4_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RK4,sys);
	
	auto success_code = predictor->Predict(RK4_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RK4_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RK4_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RK4_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
}







BOOST_AUTO_TEST_CASE(monodromy_RK4_d)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(4.641588833612776e-1), complex_dbl(7.416198487095662e-1);
	
	// Starting time
	complex_dbl current_time(0.7);
	// Time step
	complex_dbl delta_t(-0.01);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + real_mp("0.5")) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(0.412127272215145744043367969788438),
	complex_dbl(0.731436941908745436117221142349986);
	
	Vec<complex_dbl> RK4_prediction_result;
	[[maybe_unused]] double next_time;
	
	double tracking_tolerance(1e-5);
	bertini::tracking::StepMetadata meta;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RK4,sys);
	
	auto success_code = predictor->Predict(RK4_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RK4_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RK4_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RK4_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	
	
}




BOOST_AUTO_TEST_CASE(monodromy_RK4_mp)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + real_mp("0.5")) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412127272215145744043367969788438"),
	mpfr("0.731436941908745436117221142349986");
	
	Vec<mpfr> RK4_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	bertini::tracking::StepMetadata meta;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RK4,sys);
	
	auto success_code = predictor->Predict(RK4_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RK4_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RK4_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RK4_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
}























//////////////////////////////////////////////
//
//	RKF45
//
////////////////////////


BOOST_AUTO_TEST_CASE(circle_line_RKF45_double)
{
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(2.3,0.2), complex_dbl(1.1, 1.87);
	
	// Starting time
	complex_dbl current_time(0.9);
	// Time step
	complex_dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(2.39189497719010446148169962134860,0.215706089331670902152009632918759),
	complex_dbl(0.524023229576057910435628847490594, 1.42873163348439955724728985152555);
	double predicted_error = 0.0000106466724075688025735071053994891;
	
	Vec<complex_dbl> RKF45_prediction_result;

	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKF45,sys);
	
	auto success_code = predictor->Predict(RKF45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKF45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKF45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKF45_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
	
	
	
}










BOOST_AUTO_TEST_CASE(circle_line_RKF45_mp)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39189497719010446148169962134860","0.215706089331670902152009632918759"),
	mpfr("0.524023229576057910435628847490594", "1.42873163348439955724728985152555");
	auto predicted_error = NumFromString<NumErrorT>("0.0000106466724075688025735071053994891");
	
	Vec<mpfr> RKF45_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKF45,sys);
	
	auto success_code = predictor->Predict(RKF45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKF45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKF45_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKF45_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}







BOOST_AUTO_TEST_CASE(monodromy_RKF45_d)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(0.464158883361277585510862309093), complex_dbl(0.74161984870956629487113974408);
	
	// Starting time
	complex_dbl current_time(0.7);
	// Time step
	complex_dbl delta_t(-0.01);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	bertini::tracking::StepMetadata meta;
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(0.412128542780464095503570026382729),
	complex_dbl(0.731436941916416473300161135742533);
	double predicted_error = 7.17724133646795598396247354053062e-8;
	
	Vec<complex_dbl> RKF45_prediction_result;
	[[maybe_unused]] double next_time;
	
	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKF45,sys);
	
	auto success_code = predictor->Predict(RKF45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKF45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKF45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKF45_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
}



BOOST_AUTO_TEST_CASE(monodromy_RKF45_mp)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	bertini::tracking::StepMetadata meta;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412128542780464095503570026382729"),
	mpfr("0.731436941916416473300161135742533");
	auto predicted_error = NumFromString<NumErrorT>("7.17724133646795598396247354053062e-8");
	
	Vec<mpfr> RKF45_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKF45,sys);
	
	auto success_code = predictor->Predict(RKF45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKF45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKF45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKF45_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}

























//////////////////////////////////////////////
//
//	RK Cash-Karp 45
//
////////////////////////

BOOST_AUTO_TEST_CASE(circle_line_RKCK45_double)
{
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(2.3,0.2), complex_dbl(1.1, 1.87);
	
	// Starting time
	complex_dbl current_time(0.9);
	// Time step
	complex_dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(2.39189687053703334440737233377404,0.215710089694839238261207432302796),
	complex_dbl(0.524023000601737797891275060536396, 1.42873127596071439996076584113815);
	double predicted_error = 0.00000353010590253211978478006394088836;
	
	Vec<complex_dbl> RKCK45_prediction_result;

	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKCashKarp45,sys);
	
	auto success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
	
	
	
}










BOOST_AUTO_TEST_CASE(circle_line_RKCK45_mp)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39189687053703334440737233377404","0.215710089694839238261207432302796"),
	mpfr("0.524023000601737797891275060536396", "1.42873127596071439996076584113815");
	auto predicted_error = NumFromString<NumErrorT>("0.00000353010590253211978478006394088836");
	
	Vec<mpfr> RKCK45_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKCashKarp45,sys);
	
	auto success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}







BOOST_AUTO_TEST_CASE(monodromy_RKCK45_d)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(0.464158883361277585510862309093), complex_dbl(0.74161984870956629487113974408);
	
	// Starting time
	complex_dbl current_time(0.7);
	// Time step
	complex_dbl delta_t(-0.01);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	bertini::tracking::StepMetadata meta;
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(0.412128535278042242819741034030722),
	complex_dbl(0.731436941916396784391576913351911);
	double predicted_error = 4.51352044466211707817977052519894e-9;
	
	Vec<complex_dbl> RKCK45_prediction_result;
	[[maybe_unused]] double next_time;
	
	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKCashKarp45,sys);
	
	auto success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
}



BOOST_AUTO_TEST_CASE(monodromy_RKCK45_mp)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	bertini::tracking::StepMetadata meta;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412128535278042242819741034030722"),
	mpfr("0.731436941916396784391576913351911");
	auto predicted_error = NumFromString<NumErrorT>("4.51352044466211707817977052519894e-9");
	
	Vec<mpfr> RKCK45_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKCashKarp45,sys);
	
	auto success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}





BOOST_AUTO_TEST_CASE(monodromy_RKCK45_mp_change_precision)
{

	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	bertini::tracking::StepMetadata meta;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412128535278042242819741034030722"),
	mpfr("0.731436941916396784391576913351911");
	auto predicted_error = NumFromString<NumErrorT>("4.51352044466211707817977052519894e-9");
	
	Vec<mpfr> RKCK45_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKCashKarp45,sys);
	
	auto success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
	
	
	
	
	
	DefaultPrecision(50);
	
	Precision(current_space, 50);
	Precision(current_time, 50);
	Precision(delta_t, 50);

	// Starting point in spacetime step
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	current_time = mpfr("0.7");
	// Time step
	delta_t = mpfr("-0.01");
	
	
	
	
	sys.precision(50);
	AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	
	predicted << mpfr("0.41212853527804224281974103403072207383998320746093"),
	mpfr("0.73143694191639678439157691335191077020461981185497");
	predicted_error = NumFromString<NumErrorT>("4.5135204446621170781797705326691218021056435215073e-9");
	
	
	predictor->ChangePrecision(50);
	success_code = predictor->Predict(RKCK45_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKCK45_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKCK45_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKCK45_prediction_result(ii)-predicted(ii)) < 1e-47);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < 1e-47);
	
}
















//////////////////////////////////////////////
//
//	RK Dormand-Prince 56
//
////////////////////////

BOOST_AUTO_TEST_CASE(circle_line_RKDP56_double)
{
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(2.3,0.2), complex_dbl(1.1, 1.87);
	
	// Starting time
	complex_dbl current_time(0.9);
	// Time step
	complex_dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(2.39189763095027864748166494355925,0.215711752936893239277981497324557),
	complex_dbl(0.524022748677715856115185568097945, 1.42873072156957928016044855615010);
	double predicted_error = 6.79397491522542193110307157970405e-7;
	
	Vec<complex_dbl> RKDP56_prediction_result;

	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKDormandPrince56,sys);
	
	auto success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
	
	
	
}










BOOST_AUTO_TEST_CASE(circle_line_RKDP56_mp)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39189763095027864748166494355925","0.215711752936893239277981497324557"),
	mpfr("0.524022748677715856115185568097945", "1.42873072156957928016044855615010");
	auto predicted_error = NumFromString<NumErrorT>("6.79397491522542193110307157970405e-7");
	
	Vec<mpfr> RKDP56_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKDormandPrince56,sys);
	
	auto success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}








BOOST_AUTO_TEST_CASE(circle_line_RKDP56_mp_change_precision)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39189763095027864748166494355925","0.215711752936893239277981497324557"),
	mpfr("0.524022748677715856115185568097945", "1.42873072156957928016044855615010");
	auto predicted_error = NumFromString<NumErrorT>("6.79397491522542193110307157970405e-7");
	
	Vec<mpfr> RKDP56_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKDormandPrince56,sys);
	
	auto success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);

	
	
	
	
	
	
	bertini::DefaultPrecision(50);
	Precision(current_space, 50);
	Precision(current_time, 50);
	Precision(delta_t, 50);
	
	// Starting point in spacetime step
	current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
	
	// Starting time
	current_time = mpfr("0.9");
	// Time step
	delta_t = mpfr("-0.1");
	
	
	
	sys.precision(50);
	AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	
	predicted << mpfr("2.3918976309502786474816649435592524145893505795708","0.21571175293689323927798149732455717990784116340616"),
	mpfr("0.52402274867771585611518556809794390786903320453982", "1.4287307215695792801604485561500984044649241859097");
	predicted_error = NumFromString<NumErrorT>("6.7939749152254219311030715790073321381093755581241e-7");
	
	
	predictor->ChangePrecision(50);
	success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < 1e-47);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < 1e-47);

}









BOOST_AUTO_TEST_CASE(monodromy_RKDP56_d)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(0.464158883361277585510862309093), complex_dbl(0.74161984870956629487113974408);
	
	// Starting time
	complex_dbl current_time(0.7);
	// Time step
	complex_dbl delta_t(-0.01);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	bertini::tracking::StepMetadata meta;
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(0.412128532164122346459968880922735),
	complex_dbl(0.731436941916392989685864031055020);
	double predicted_error = 3.85904197101299548102733617445410e-9;
	
	Vec<complex_dbl> RKDP56_prediction_result;
	[[maybe_unused]] double next_time;
	
	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKDormandPrince56,sys);
	
	auto success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
}



BOOST_AUTO_TEST_CASE(monodromy_RKDP56_mp)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	bertini::tracking::StepMetadata meta;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412128532164122346459968880922735"),
	mpfr("0.731436941916392989685864031055020");
	auto predicted_error = NumFromString<NumErrorT>("3.85904197101299548102733617445410e-9");
	
	Vec<mpfr> RKDP56_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKDormandPrince56,sys);
	
	auto success_code = predictor->Predict(RKDP56_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKDP56_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKDP56_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKDP56_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
	
}












//////////////////////////////////////////////
//
//	RK Verner 67
//
////////////////////////

BOOST_AUTO_TEST_CASE(circle_line_RKV67_double)
{
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(2.3,0.2), complex_dbl(1.1, 1.87);
	
	// Starting time
	complex_dbl current_time(0.9);
	// Time step
	complex_dbl delta_t(-0.1);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(2.39189815934576660586899846426669,0.215712712024488132602524062094065),
	complex_dbl(0.524022631256496309806889230162953, 1.42873050843900263719943909731242);
	double predicted_error = 0.00000128891520195955347062706145253149;
	
	Vec<complex_dbl> RKV67_prediction_result;

	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKVerner67,sys);
	
	auto success_code = predictor->Predict(RKV67_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKV67_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKV67_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKV67_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
	
	
	
	
}










BOOST_AUTO_TEST_CASE(circle_line_RKV67_mp)
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
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	bertini::tracking::StepMetadata meta;
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("2.39189815934576660586899846426669","0.215712712024488132602524062094065"),
	mpfr("0.524022631256496309806889230162953", "1.42873050843900263719943909731242");
	auto predicted_error = NumFromString<NumErrorT>("0.00000128891520195955347062706145253149");
	
	Vec<mpfr> RKV67_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKVerner67,sys);
	
	auto success_code = predictor->Predict(RKV67_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKV67_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKV67_prediction_result.size(); ++ii)
		BOOST_CHECK(abs(RKV67_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
}







BOOST_AUTO_TEST_CASE(monodromy_RKV67_d)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<complex_dbl> current_space(2);
	current_space << complex_dbl(0.464158883361277585510862309093), complex_dbl(0.74161984870956629487113974408);
	
	// Starting time
	complex_dbl current_time(0.7);
	// Time step
	complex_dbl delta_t(-0.01);
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	bertini::tracking::StepMetadata meta;
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	
	Vec<complex_dbl> predicted(2);
	predicted << complex_dbl(0.412128533889452110491490000899263),
	complex_dbl(0.731436941916389669876029584806957);
	double predicted_error = 1.42794733055750714441060080061e-8;
	
	Vec<complex_dbl> RKV67_prediction_result;
	[[maybe_unused]] double next_time;
	
	double tracking_tolerance(1e-5);
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKVerner67,sys);
	
	auto success_code = predictor->Predict(RKV67_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKV67_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKV67_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKV67_prediction_result(ii)-predicted(ii)) < threshold_clearance_d);
	}
	
	BOOST_CHECK(fabs(meta.error_estimate - predicted_error) < threshold_clearance_d);
}



BOOST_AUTO_TEST_CASE(monodromy_RKV67_mp)
{
	DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
	
	// Starting point in spacetime step
	Vec<mpfr> current_space(2);
	current_space << mpfr("0.464158883361277585510862309093"), mpfr("0.74161984870956629487113974408");
	
	// Starting time
	mpfr current_time("0.7");
	// Time step
	mpfr delta_t("-0.01");
	
	
	
	
	bertini::System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
	std::shared_ptr<Complex> half = Complex::Make("0.5");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,3)-1) + (1-t)*(pow(x,3) + 2) );
	sys.AddFunction( t*(pow(y,2)-1) + (1-t)*(pow(y,2) + half) );
	
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,3);
	AMP.coefficient_bound = 2;
	
	bertini::tracking::StepMetadata meta;
	
	
	Vec<mpfr> predicted(2);
	predicted << mpfr("0.412128533889452110491490000899263"),
	mpfr("0.731436941916389669876029584806957");
	auto predicted_error = NumFromString<NumErrorT>("1.42794733055750714441060080061e-8");
	
	Vec<mpfr> RKV67_prediction_result;
	mpfr next_time;
	
	double tracking_tolerance = 1e-5;
	unsigned num_steps_since_last_condition_number_computation = 1;
	unsigned frequency_of_CN_estimation = 1;
	
	std::shared_ptr<ExplicitRKPredictor> predictor = std::make_shared< ExplicitRKPredictor >(bertini::tracking::Predictor::RKVerner67,sys);
	
	auto success_code = predictor->Predict(RKV67_prediction_result, meta, sys, current_space, current_time, delta_t, num_steps_since_last_condition_number_computation, frequency_of_CN_estimation, tracking_tolerance, &AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(RKV67_prediction_result.size(),2);
	for (unsigned ii = 0; ii < RKV67_prediction_result.size(); ++ii)
	{
		BOOST_CHECK(abs(RKV67_prediction_result(ii)-predicted(ii)) < threshold_clearance_mp);
	}
	
	BOOST_CHECK(abs(meta.error_estimate - predicted_error) < threshold_clearance_mp);
}


// --- size_proportion / error_estimate behaviour (AMP precision-escalation work) -------------------
// The default predictor was Euler (no error estimate); its size_proportion fallback
// maxCoeff(K)/|delta_t|^p blows up as the step shrinks (~1/|delta_t|), spuriously inflating AMP
// precision (DigitsB).  RKF45 has an embedded error estimate, so size_proportion =
// meta.error_estimate/|delta_t|^(p+1) is a bounded, stable "a" coefficient.  These tests pin that and guard the
// default.  (Previously size_proportion / error_estimate were computed in tests but never asserted.)

namespace {
	// build the smooth test homotopy + a generic on-path point used by the size_proportion tests.
	struct PredFixture {
		bertini::System sys;
		Vec<complex_dbl> current_space{2};
		complex_dbl current_time{0.9};
		bertini::tracking::AdaptiveMultiplePrecisionConfig AMP;
		PredFixture() {
			Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");
			sys.AddVariableGroup(VariableGroup{x,y});
			sys.AddPathVariable(t);
			sys.AddFunction( t*(pow(x,2)-1) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
			sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
			current_space << complex_dbl(2.3,0.2), complex_dbl(1.1,1.87);
			AMP = bertini::tracking::AMPConfigFrom(sys);
			AMP.coefficient_bound = 5;
		}
	};
}

BOOST_AUTO_TEST_CASE(default_predictor_is_rkf45)
{
	// guards the decision to default to an error-estimate predictor (was Euler, atrocious for
	// performance and the cause of the size_proportion blow-up).
	BOOST_CHECK(bertini::tracking::predict::DefaultPredictor() == bertini::tracking::Predictor::RKF45);
}

BOOST_AUTO_TEST_CASE(rkf45_size_proportion_stays_bounded_as_step_shrinks)
{
	PredFixture f;
	auto predictor = std::make_shared<ExplicitRKPredictor>(bertini::tracking::Predictor::RKF45, f.sys);

	double tracking_tolerance(1e-5); unsigned nsc(1), freq(1);
	Vec<complex_dbl> result; bertini::tracking::StepMetadata meta;

	double sp_max = 0.0, sp_min = 1e300;
	for (double h : {-0.1, -0.05, -0.025, -0.0125}) // step shrinks 8x; meta.error_estimate stays above roundoff
	{
		auto code = predictor->Predict(result, meta, f.sys, f.current_space, f.current_time, complex_dbl(h), nsc, freq, tracking_tolerance, &f.AMP);
		BOOST_REQUIRE(code == bertini::SuccessCode::Success);
		BOOST_REQUIRE(meta.size_proportion > 0.0);
		sp_max = std::max(sp_max, meta.size_proportion);
		sp_min = std::min(sp_min, meta.size_proportion);
	}
	// RKF45's size_proportion ~ constant: it stays within a small factor as the step shrinks 8x.
	// Euler's fallback maxCoeff(K)/|delta_t| would vary by ~8x over the same range (and unboundedly
	// as the step keeps shrinking) -- which is what spuriously escalated precision.
	BOOST_CHECK_LT(sp_max / sp_min, 5.0);
}

BOOST_AUTO_TEST_CASE(rkf45_error_estimate_has_order_p_plus_1)
{
	PredFixture f;
	auto predictor = std::make_shared<ExplicitRKPredictor>(bertini::tracking::Predictor::RKF45, f.sys);
	const unsigned p = bertini::tracking::predict::Order(bertini::tracking::Predictor::RKF45); // 4

	double tracking_tolerance(1e-5); unsigned nsc(1), freq(1);
	Vec<complex_dbl> result; bertini::tracking::StepMetadata meta;

	auto err_at = [&](double h) {
		auto code = predictor->Predict(result, meta, f.sys, f.current_space, f.current_time, complex_dbl(h), nsc, freq, tracking_tolerance, &f.AMP);
		BOOST_REQUIRE(code == bertini::SuccessCode::Success);
		return meta.error_estimate;
	};

	// error_estimate ~ C * h^(p+1): halving the step divides the estimate by 2^(p+1).
	const double expected_ratio = std::pow(2.0, double(p + 1)); // 2^5 = 32 for RKF45
	for (double h : {-0.1, -0.05, -0.025}) // above the roundoff floor in double
	{
		double ratio = err_at(h) / err_at(h / 2.0);
		// within 2x of the asymptotic ratio (higher-order terms perturb it a little)
		BOOST_CHECK_GT(ratio, expected_ratio / 2.0);
		BOOST_CHECK_LT(ratio, expected_ratio * 2.0);
	}
}


BOOST_AUTO_TEST_CASE(euler_size_proportion_stays_bounded_as_step_shrinks)
{
	// Regression for the non-error-estimate (Euler) size_proportion bug: AMP2 (Eqs 9-10) gives
	// a = ||d||/|s| = ||K|| ~ maxCoeff(|K|), an O(1) step-independent constant.  The old fallback
	// maxCoeff(K)/|delta_t|^p over-divided by the step, so $a$ grew ~1/|delta_t| as the step shrank
	// and spuriously escalated AMP precision.  The corrected formula must stay (nearly) constant.
	PredFixture f;
	auto predictor = std::make_shared<ExplicitRKPredictor>(bertini::tracking::Predictor::Euler, f.sys);

	double tracking_tolerance(1e-5); unsigned nsc(1), freq(1);
	Vec<complex_dbl> result; bertini::tracking::StepMetadata meta;

	double sp_max = 0.0, sp_min = 1e300;
	for (double h : {-0.1, -0.05, -0.025, -0.0125}) // step shrinks 8x
	{
		// Euler has no error estimate, so use the size_proportion-only Predict overload.
		auto code = predictor->Predict(result, meta, f.sys, f.current_space, f.current_time, complex_dbl(h), nsc, freq, tracking_tolerance, &f.AMP);
		BOOST_REQUIRE(code == bertini::SuccessCode::Success);
		BOOST_REQUIRE(meta.size_proportion > 0.0);
		sp_max = std::max(sp_max, meta.size_proportion);
		sp_min = std::min(sp_min, meta.size_proportion);
	}
	// ||K|| (the tangent magnitude) barely moves over a smooth path as the step shrinks 8x; the
	// buggy maxCoeff(K)/|delta_t| would vary by ~8x over the same range (and unboundedly as h->0).
	BOOST_CHECK_LT(sp_max / sp_min, 2.0);
}

BOOST_AUTO_TEST_SUITE_END()
















