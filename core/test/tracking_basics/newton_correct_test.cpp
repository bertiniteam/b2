//This file is part of Bertini 2.
//
//newton_correct_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//newton_correct_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with newton_correct_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
#include "trackers/newton_corrector.hpp"




extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;



BOOST_AUTO_TEST_SUITE(newton_correct_tracking_basics)

using System = bertini::System;
using Variable = bertini::node::Variable;
using bertini::MakeVariable;
using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;
using NewtonCorrector = bertini::tracking::correct::NewtonCorrector;



using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;
using mpq_rational = bertini::mpq_rational;

template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

using bertini::DefaultPrecision;



BOOST_AUTO_TEST_CASE(circle_line_one_corrector_step_double)
{
	
	// Starting point in spacetime step
	Vec<dbl> current_space(2);
	current_space << dbl(2.3,0.2), dbl(1.1, 1.87);
	
	// Starting time
	dbl current_time(0.9);
	// Time step
	
	
	bertini::System sys;
	Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
	
	VariableGroup vars{x,y};
	
	sys.AddVariableGroup(vars);
	sys.AddPathVariable(t);
	
	// Define homotopy system
	sys.AddFunction( t*(pow(x,2)-1.0) + (1-t)*(pow(x,2) + pow(y,2) - 4) );
	sys.AddFunction( t*(y-1) + (1-t)*(2*x + 5*y) );
	
	
	auto AMP = bertini::tracking::AMPConfigFrom(sys);
	
	BOOST_CHECK_EQUAL(AMP.degree_bound,2);
	AMP.coefficient_bound = 5;
	
	
	
	double tracking_tolerance(1e-5);
	
	
	Vec<dbl> corrected(2);
	corrected << dbl(1.36296628875178620892887063382866,0.135404746200380445814213878747082),
	dbl(0.448147673035459113010161338024478, -0.0193435351714829208306019826781546);
	
	Vec<dbl> newton_correction_result;
	
	tracking_tolerance = 1e1;
	unsigned max_num_newton_iterations = 1;
	unsigned min_num_newton_iterations = 1;
	std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);

	auto success_code = corrector->Correct(newton_correction_result,
											  sys,
											  current_space,
											  current_time,
											  tracking_tolerance,
											  min_num_newton_iterations,
											  max_num_newton_iterations,
											  AMP);
	
	BOOST_CHECK(success_code==bertini::SuccessCode::Success);
	BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
	for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_d);
	
	}
	
	BOOST_AUTO_TEST_CASE(circle_line_one_corrector_step_mp)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
		
		// Starting time
		mpfr current_time("0.9");
		// Time step
		
		
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
		
		
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("1.36296628875178620892887063382866","0.135404746200380445814213878747082"),
		mpfr("0.448147673035459113010161338024478", "-0.0193435351714829208306019826781546");
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
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
		
		
		
		double tracking_tolerance(1e-5);
		
		
		Vec<dbl> corrected(2);
		corrected << dbl(1.14542104415948767661671388923986, 0.0217584797792294631577151109622764),
		dbl(0.47922556512007318905475515868002, -0.00310835425417563759395930156603948);
		
		Vec<dbl> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 2;
		unsigned min_num_newton_iterations = 2;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
		{
			BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_d);
		}
		
		
	}
	
	BOOST_AUTO_TEST_CASE(circle_line_two_corrector_steps_mp)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		// Starting point in spacetime step
		Vec<mpfr> current_space(2);
		current_space << mpfr("2.3","0.2"), mpfr("1.1", "1.87");
		
		// Starting time
		mpfr current_time("0.9");
		// Time step
		
		
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
		
		
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("1.14542104415948767661671388923986", "0.0217584797792294631577151109622764"),
		mpfr("0.47922556512007318905475515868002", "-0.00310835425417563759395930156603948");
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 2;
		unsigned min_num_newton_iterations = 2;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK(success_code==bertini::SuccessCode::Success);
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		for (unsigned ii = 0; ii < newton_correction_result.size(); ++ii)
			BOOST_CHECK(abs(newton_correction_result(ii)-corrected(ii)) < threshold_clearance_mp);
		
	}
	
	BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_B_violated_double)
	{
		
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence.
		 Have trimmed current_space to have 16 digits after decimal point. Also, saftey_digits_1 has been set to 32000
		 to set off the AMPCriterionB condition.
		 */
		Vec<dbl> current_space(2);
		current_space << dbl(25.4088532364492429,-38.0519122331723744),
		dbl(-0.0212298348984663,-0.1778146465316983);
		
		dbl current_time(0);
		dbl delta_t(.1);
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 5;
		AMP.safety_digits_1 = 32000;
		
		
		Vec<dbl> newton_correction_result;
		
		double tracking_tolerance = double(10);
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_B_violated_mp)
	{
		
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence.
		 Also, saftey_digits_1 has been set to 32000 to set off the AMPCriterionB condition.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("25.408853236449242927412","-38.051912233172374487976"),
		mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");
		
		mpfr current_time("0");
		mpfr delta_t(".1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 5;
		AMP.safety_digits_1 = 32000;
		
		double tracking_tolerance = 1e-5;
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_C_violated_double)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence.
		 Have trimmed current_space to have 16 digits after decimal point. Also, saftey_digits_2 has been set to 32000
		 to set off the AMPCriterionC condition.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("256185069753.4088532364492429","-387520022558.0519122331723744"),
		mpfr("-0.0212298348984663","-0.1778146465316983");
		
		mpfr current_time("0");
		mpfr delta_t(".1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 5;
		AMP.safety_digits_2 = 32000;
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("3701884101067.778","-5599679215240.413"),
		mpfr("-1.043206463433583e25","-2.450083921191992e25");
		
		
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_amp_criterion_C_violated_mp)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence.
		 Also, saftey_digits_2 has been set to 32000 to set off the AMPCriterionC condition.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
		mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");
		
		mpfr current_time("0");
		mpfr delta_t(".1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,3);
		AMP.coefficient_bound = 5;
		AMP.safety_digits_2 = 32000;
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("3701884101067.778","-5599679215240.413"),
		mpfr("-1.043206463433583e25","-2.450083921191992e25");
		
		
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::HigherPrecisionNecessary);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_linear_algebra_fails_double)
	{
		
		/*Test case checking to make sure a polynomial system that is 0 will fail when the newton correction step is done.
		 This test case differs from above as it has the precision type set to double.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("0","0"),mpfr("0","0");
		
		mpfr current_time("1");
		mpfr delta_t("-.1");
		
		
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpfr("0","0")*x);
		sys.AddFunction(mpfr("0","0")*y);
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,1);
		AMP.coefficient_bound = 5;
		
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("0","0"),mpfr("0","0");
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::MatrixSolveFailure);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_linear_algebra_fails_mp)
	{
		
		/*Test case checking to make sure a polynomial system that is 0 will fail when the newton correction step is done.
		 This test case differs from above as it has the precision type set to adaptive.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("0","0"),mpfr("0","0");
		
		mpfr current_time("1");
		mpfr delta_t("-.1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpfr("0","0")*x);
		sys.AddFunction(mpfr("0","0")*y);
		
		auto AMP = bertini::tracking::AMPConfigFrom(sys);
		
		BOOST_CHECK_EQUAL(AMP.degree_bound,1);
		AMP.coefficient_bound = 5;
		
		
		double tracking_tolerance = 1e-5;
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("0","0"),mpfr("0","0");
		
		Vec<mpfr> newton_correction_result;
		
		tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations,
												  AMP);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::MatrixSolveFailure);
	}
	
	
	BOOST_AUTO_TEST_CASE(newton_step_diverging_to_infinity_fails_to_converge_d)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
		mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");
		
		mpfr current_time("0");
		mpfr delta_t(".1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("3701884101067.778","-5599679215240.413"),
		mpfr("-1.043206463433583e25","-2.450083921191992e25");
		
		
		double tracking_tolerance = 1e1;
		
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		
		Vec<mpfr> newton_correction_result;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::FailedToConverge);
	}
	
	BOOST_AUTO_TEST_CASE(newton_step_diverging_to_infinity_fails_to_converge_mp)
	{
		DefaultPrecision(TRACKING_TEST_MPFR_DEFAULT_DIGITS);
		/*
		 Using the Griewank Osborne example. Starting at t = 0 where there is a multiplicity 3 isolated solution. We predict
		 to .1 and try to correct back down. Anywhere except at t = 0, we will have divergence. The difference from this
		 and the above version is we will use PrecisionType = Adaptive.
		 */
		Vec<mpfr> current_space(2);
		current_space << mpfr("256185069753.408853236449242927412","-387520022558.051912233172374487976"),
		mpfr("-0.0212298348984663761753389403711889","-0.177814646531698303094367623155171");
		
		mpfr current_time("0");
		mpfr delta_t(".1");
		
		current_time += delta_t;
		
		bertini::System sys;
		Var x = MakeVariable("x"), y = MakeVariable("y"), t = MakeVariable("t");
		
		VariableGroup vars{x,y};
		
		sys.AddVariableGroup(vars);
		sys.AddPathVariable(t);
		
		sys.AddFunction(mpq_rational(29,16)*pow(x,3) - 2*x*y + t);
		sys.AddFunction(y - pow(x,2));
		
		
		Vec<mpfr> corrected(2);
		corrected << mpfr("3701884101067.778","-5599679215240.413"),
		mpfr("-1.043206463433583e25","-2.450083921191992e25");
		
		
		
		Vec<mpfr> newton_correction_result;
		
		double tracking_tolerance = 1e1;
		unsigned max_num_newton_iterations = 1;
		unsigned min_num_newton_iterations = 1;
		std::shared_ptr<NewtonCorrector> corrector = std::make_shared<NewtonCorrector>(sys);
		auto success_code = corrector->Correct(newton_correction_result,
												  sys,
												  current_space,
												  current_time,
												  tracking_tolerance,
												  min_num_newton_iterations,
												  max_num_newton_iterations);
		
		BOOST_CHECK_EQUAL(newton_correction_result.size(),2);
		BOOST_CHECK(success_code==bertini::SuccessCode::FailedToConverge);
	}

BOOST_AUTO_TEST_SUITE_END()





