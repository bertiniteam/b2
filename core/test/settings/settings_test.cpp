//This file is part of Bertini 2.
//
//test/settings/settings_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/settings/settings_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/settings/settings_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

//  Created by Collins, James B. on 8/27/15.
//  Copyright (c) 2015 West Texas A&M University.
//
// additionally authored by Dani Brake, University of Notre Dame

//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Settings Testing"


#include <iostream>

#include <cstdlib>
#include <cmath>

#include "bertini2/bertini.hpp"
#include "bertini2/function_tree.hpp"
#include <bertini2/io/parsing/settings_parsers.hpp>
#include <bertini2/io/parsing/classic_utilities.hpp>


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

using mpfr = bertini::mpfr;
using mpfr_float = bertini::mpfr_float;
using mpq_rational = bertini::mpq_rational;

namespace algorithm = bertini::algorithm;

using std::abs;

BOOST_AUTO_TEST_SUITE(config_settings)

BOOST_AUTO_TEST_CASE(read_mptype)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n MPType  : 0; \n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	PrecisionType type;
	ConfigSettingParser<std::string::const_iterator, PrecisionType> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, type);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(type == PrecisionType::Fixed);
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n  \n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, type);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(type == PrecisionType::Adaptive);
}



BOOST_AUTO_TEST_CASE(read_predictor)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n MPType: 0; \n  ODEPredictor: 4; % the predictor type\n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	Predictor predictor;
	ConfigSettingParser<std::string::const_iterator, Predictor> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == Predictor::RKNorsett34);
	
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n MPType: 0; \n  ODEPredictor : -1; end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == Predictor::Constant);
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n MPType: 0;  end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == Predictor::RKF45);
	
	
}



BOOST_AUTO_TEST_CASE(read_security)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	
	double tol = 1e-16;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n   SecurityLevel: 1;\n SecurityMaxNorm: -2.34; % the predictor type\n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	using SecurityConfig = bertini::endgame::SecurityConfig;

	SecurityConfig settings;
	ConfigSettingParser<std::string::const_iterator, SecurityConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.level == 1);
	BOOST_CHECK(abs(settings.max_norm - (-2.34)) < tol);
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.level == 1);
	BOOST_CHECK(abs(settings.max_norm - 2.34) < tol);
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n ;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	settings = SecurityConfig();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.level == 0);
	BOOST_CHECK(abs(settings.max_norm - 2.34) < tol);
	
	
	
	inputfile = ParseInputFile("Config \n heLlo: 9 \n  SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	settings = SecurityConfig();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.level == 1);
	BOOST_CHECK(settings.max_norm == 10000);
	
	inputfile = ParseInputFile("Config \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	settings = SecurityConfig();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.level == 0);
	BOOST_CHECK(settings.max_norm == 10000);
	
}





BOOST_AUTO_TEST_CASE(read_tolerances)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	bertini::mpfr_float::default_precision(30);
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n   FinalTol   : -.845e-7   ;\n TrackTolDuringEG   : 234e-4   ; % the predictor type\n TrackTolBeforeEG         : 7.32e3; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	algorithm::TolerancesConfig tols;
	ConfigSettingParser<std::string::const_iterator, algorithm::TolerancesConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, tols);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(abs(tols.newton_before_endgame -  7.32e3) < tol);
	BOOST_CHECK(abs(tols.newton_during_endgame -  234e-4) < tol);
	BOOST_CHECK(abs(tols.final_tolerance -  -0.845e-7) < tol);
	BOOST_CHECK(abs(tols.path_truncation_threshold -  100000) < tol);
	
	
}



BOOST_AUTO_TEST_CASE(read_stepping)
{
	bool check_rational = false;

	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n StepSuccessFactor  : 4.2;  FinalTol: 1.845e-7;\n MaxNumberSteps: 234; % the predictor type\nMaxStepSize: 1e-2; StepsForIncrease: 7;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	SteppingConfig structure;
	ConfigSettingParser<std::string::const_iterator,SteppingConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);

	if (check_rational)
	{
		BOOST_CHECK_EQUAL(structure.max_step_size, mpq_rational(1,100));
		BOOST_CHECK_EQUAL(structure.step_size_success_factor, mpq_rational(42,10));
		BOOST_CHECK_EQUAL(structure.step_size_fail_factor, mpq_rational(1/2));
	}
	else
	{
		double tol_double = 1e-14;
		BOOST_CHECK_CLOSE(structure.max_step_size, mpq_rational(1,100), tol_double);
		BOOST_CHECK_CLOSE(structure.step_size_success_factor, mpq_rational(42,10), tol_double);
		BOOST_CHECK_CLOSE(structure.step_size_fail_factor, mpq_rational(1/2), tol_double);
	}

	BOOST_CHECK_EQUAL(structure.consecutive_successful_steps_before_stepsize_increase, 7);
	BOOST_CHECK_EQUAL(structure.max_num_steps, 234);
	
	
}





BOOST_AUTO_TEST_CASE(read_newton)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n MaxNewtonIts  : 5; \n heLlo: 9 \n StepSuccessFactor: 4.2;   end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	NewtonConfig structure;
	ConfigSettingParser<std::string::const_iterator,NewtonConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_num_newton_iterations == 5);
	
	
	inputfile = ParseInputFile("Config \n  \n heLlo: 9 \n StepSuccessFactor: 4.2;  \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	structure = NewtonConfig();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_num_newton_iterations == 2);
	
}






BOOST_AUTO_TEST_CASE(read_endgame)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	bertini::mpfr_float::default_precision(30);
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n   NumSamplePoints: 34;\n TrackTolDuringEG: 234e-4; % the predictor type\n NbhdRadius: 4.3e-7; \n SampleFactor   : 8e-3;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	bertini::endgame::EndgameConfig settings;
	ConfigSettingParser<std::string::const_iterator,bertini::endgame::EndgameConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(settings.num_sample_points == 34);
	BOOST_CHECK(abs(settings.min_track_time - 4.3e-7) < tol);
	BOOST_CHECK_CLOSE(settings.sample_factor,  mpq_rational(8,1000), tol);
	
	
}





BOOST_AUTO_TEST_CASE(read_powerseries)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n MaxNewtonIts: 5; \n heLlo: 9 \n MaxCycleNum : 4; StepSuccessFactor: 4.2;   end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	bertini::endgame::PowerSeriesConfig structure;
	ConfigSettingParser<std::string::const_iterator,bertini::endgame::PowerSeriesConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_cycle_number == 4);
	
	
	inputfile = ParseInputFile("Config \n  \n heLlo: 9 \n StepSuccessFactor: 4.2;  \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	structure = bertini::endgame::PowerSeriesConfig();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_cycle_number == 6);
	
}




BOOST_AUTO_TEST_CASE(read_cauchy)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = ParseInputFile("Config \n heLlo: 9 \n StepSuccessFactor: 4.2;  CycleTimeCutoff: 5.76e2  ;\n MaxNumberSteps: 234; % the predictor type \n RAtioTimeCutoff: 1e-12; StepsForIncrease: 7;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	bertini::endgame::CauchyConfig structure;
	ConfigSettingParser<std::string::const_iterator,bertini::endgame::CauchyConfig> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(abs(structure.cycle_cutoff_time - 5.76e2) < tol);
	BOOST_CHECK(abs(structure.ratio_cutoff_time / 1e-12 - 1) < tol);
	
	
}



BOOST_AUTO_TEST_CASE(read_meta)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::algorithm;

	
	
	SplitInputFile inputfile = ParseInputFile("Config \n TrackType: 5; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	bertini::algorithm::classic::AlgoChoice tt;
	ConfigSettingParser<std::string::const_iterator,bertini::algorithm::classic::AlgoChoice> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, tt);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(tt == bertini::algorithm::classic::AlgoChoice::WitnessSetProjection);
	
	
	inputfile = ParseInputFile("Config \n  \n heLlo: 9 \n StepSuccessFactor: 4.2;  \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	tt = bertini::algorithm::classic::AlgoChoice();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, tt);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(tt == bertini::algorithm::classic::AlgoChoice::ZeroDim);
	
}



BOOST_AUTO_TEST_CASE(all_config_settings)
{
	using namespace bertini::parsing::classic;
	using namespace bertini::tracking;
	double tol = 1e-15;
	
	SplitInputFile inputfile = ParseInputFile("Config \n ODEPredictor: 7; \n MPType: 0; \n MaxNewtonIts: 7;  FinalTol: 1.845e-7;\n SampleFactor: 0.647; \n NumSamplePoints: 7;\n\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	auto sets = ConfigParser<Predictor, NewtonConfig, algorithm::TolerancesConfig, bertini::endgame::EndgameConfig, SteppingConfig>::Parse(configStr);
	
	SteppingConfig steps = std::get<SteppingConfig>(sets);
	Predictor pred = std::get<Predictor>(sets);
	NewtonConfig newt = std::get<NewtonConfig>(sets);
	algorithm::TolerancesConfig tols = std::get<algorithm::TolerancesConfig>(sets);
	bertini::endgame::EndgameConfig end = std::get<bertini::endgame::EndgameConfig>(sets);
	
	BOOST_CHECK(pred == Predictor::RKDormandPrince56);
	BOOST_CHECK_EQUAL( steps.max_step_size, mpq_rational(1,10));
	BOOST_CHECK_EQUAL(newt.max_num_newton_iterations, 7);
	BOOST_CHECK_EQUAL(newt.min_num_newton_iterations, 1);
	BOOST_CHECK(abs(tols.final_tolerance - 1.845e-7) < tol);
	BOOST_CHECK_CLOSE( end.sample_factor, mpq_rational(647,1000), tol);
	BOOST_CHECK_EQUAL(end.num_sample_points, 7);
	BOOST_CHECK_EQUAL(end.min_track_time, 1e-100);
}



BOOST_AUTO_TEST_SUITE_END()
