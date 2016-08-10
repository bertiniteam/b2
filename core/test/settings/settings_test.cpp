//This file is part of Bertini 2.
//
//settings_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//settings_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with settings_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

//  Created by Collins, James B. on 8/27/15.
//  Copyright (c) 2015 West Texas A&M University.
//

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
#include <bertini2/classic/split_parsing.hpp>
#include <bertini2/classic/config_parsing.hpp>


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

using mpfr = bertini::mpfr;
using mpfr_float = bertini::mpfr_float;




BOOST_AUTO_TEST_SUITE(config_settings)

BOOST_AUTO_TEST_CASE(read_mptype)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n MPType: 0; \n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	PrecisionType type;
	parsing::ConfigSettingParser<std::string::const_iterator, PrecisionType> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, type);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(type == PrecisionType::Fixed);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  \n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, type);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(type == PrecisionType::Adaptive);
}



BOOST_AUTO_TEST_CASE(read_predictor)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n MPType: 0; \n  ODEPredictor: 4; % the predictor type\n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Predictor predictor;
	parsing::ConfigSettingParser<std::string::const_iterator, config::Predictor> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == config::Predictor::RKNorsett34);
	
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n MPType: 0; \n  ODEPredictor: -1; end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == config::Predictor::Constant);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n MPType: 0;  end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	
	
	iter = configStr.begin();
	end = configStr.end();
	
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, predictor);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(predictor == config::Predictor::RKF45);
	
	
}



BOOST_AUTO_TEST_CASE(read_security_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = double;
	
	double tol = 1e-16;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n   SecurityLevel: 1;\n SecurityMaxNorm: 2.34; % the predictor type\n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Security<T> security;
	parsing::ConfigSettingParser<std::string::const_iterator, config::Security<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(fabs(security.max_norm - 2.34) < tol);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(fabs(security.max_norm - 2.34) < tol);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n ;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	security = config::Security<double>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 0);
	BOOST_CHECK(fabs(security.max_norm - 2.34) < tol);
	
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	security = config::Security<double>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(security.max_norm == 100000);
	
	inputfile = parsing::ParseInputFile("Config \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	security = config::Security<double>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 0);
	BOOST_CHECK(security.max_norm == 100000);
	
}





BOOST_AUTO_TEST_CASE(read_security_mp)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	bertini::mpfr_float::default_precision(30);
	
	using T = bertini::mpfr_float;
	
	T tol = 1e-27;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n   SecurityLevel: 1;\n SecurityMaxNorm: 2.34; % the predictor type\n NeWTon: 1; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Security<T> security;
	parsing::ConfigSettingParser<std::string::const_iterator, config::Security<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(fabs(security.max_norm - mpfr_float("2.34")) < tol);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(fabs(security.max_norm - mpfr_float("2.34")) < tol);
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityMaxNorm: 2.34;\n ;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	
	security = config::Security<T>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 0);
	BOOST_CHECK(fabs(security.max_norm - mpfr_float("2.34")) < tol);
	
	
	
	inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n  SecurityLevel: 1;end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	security = config::Security<T>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 1);
	BOOST_CHECK(fabs(security.max_norm - mpfr_float("100000")) < tol);
	
	inputfile = parsing::ParseInputFile("Config \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	security = config::Security<T>();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.level == 0);
	BOOST_CHECK(abs(security.max_norm - mpfr_float("100000")) < tol);
	
}




BOOST_AUTO_TEST_CASE(read_tolerance_mp)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = mpfr_float;
	bertini::mpfr_float::default_precision(30);
	
	
	double tol = 1e-27;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n   FinalTol: 1.845e-7;\n TrackTolDuringEG: 234e-4; % the predictor type\n TrackTolBeforeEG: 7.32e3; \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Tolerances<T> security;
	parsing::ConfigSettingParser<std::string::const_iterator, config::Tolerances<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(abs(security.newton_before_endgame - mpfr_float("7.32e3")) < tol);
	BOOST_CHECK(abs(security.newton_during_endgame - mpfr_float("234e-4")) < tol);
	BOOST_CHECK(abs(security.final_tolerance - mpfr_float("1.845e-7")) < tol);
	BOOST_CHECK(abs(security.path_truncation_threshold - mpfr_float("100000")) < tol);
	
	
}



BOOST_AUTO_TEST_CASE(read_stepping_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n StepSuccessFactor: 4.2;  FinalTol: 1.845e-7;\n MaxNumberSteps: 234; % the predictor type\nMaxStepSize: 1e-2; StepsForIncrease: 7;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Stepping<T> structure;
	parsing::ConfigSettingParser<std::string::const_iterator,config::Stepping<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(fabs(structure.max_step_size - 1e-2) < tol);
	BOOST_CHECK(fabs(structure.step_size_success_factor - 4.2) < tol);
	BOOST_CHECK(fabs(structure.step_size_fail_factor - 0.5) < tol);
	BOOST_CHECK(structure.consecutive_successful_steps_before_stepsize_increase == 7);
	BOOST_CHECK(structure.max_num_steps == 234);
	
	
}





BOOST_AUTO_TEST_CASE(read_newton_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n MaxNewtonIts: 5; \n heLlo: 9 \n StepSuccessFactor: 4.2;   end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Newton structure;
	parsing::ConfigSettingParser<std::string::const_iterator,config::Newton> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_num_newton_iterations == 5);
	
	
	inputfile = parsing::ParseInputFile("Config \n  \n heLlo: 9 \n StepSuccessFactor: 4.2;  \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	structure = config::Newton();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_num_newton_iterations == 2);
	
}






BOOST_AUTO_TEST_CASE(read_endgame_mp)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = mpfr_float;
	bertini::mpfr_float::default_precision(30);
	
	
	double tol = 1e-27;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n   NumSamplePoints: 34;\n TrackTolDuringEG: 234e-4; % the predictor type\n NbhdRadius: 4.3e-7; \n SampleFactor: 8e-3;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Endgame<T> security;
	parsing::ConfigSettingParser<std::string::const_iterator,config::Endgame<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, security);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(security.num_sample_points == 34);
	BOOST_CHECK(abs(security.min_track_time - mpfr_float("4.3e-7")) < tol);
	BOOST_CHECK(abs(security.sample_factor - mpfr_float("8e-3")) < tol);
	
	
}





BOOST_AUTO_TEST_CASE(read_powerseries_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n MaxNewtonIts: 5; \n heLlo: 9 \n MaxCycleNum: 4; StepSuccessFactor: 4.2;   end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::PowerSeries structure;
	parsing::ConfigSettingParser<std::string::const_iterator,config::PowerSeries> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_cycle_number == 4);
	
	
	inputfile = parsing::ParseInputFile("Config \n  \n heLlo: 9 \n StepSuccessFactor: 4.2;  \n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	configStr = inputfile.Config();
	
	iter = configStr.begin();
	end = configStr.end();
	structure = config::PowerSeries();
	parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(structure.max_cycle_number == 6);
	
}




BOOST_AUTO_TEST_CASE(read_cauchy_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	using T = double;
	
	
	double tol = 1e-15;
	SplitInputFile inputfile = parsing::ParseInputFile("Config \n heLlo: 9 \n StepSuccessFactor: 4.2;  CycleTimeCutoff: 5.76e2;\n MaxNumberSteps: 234; % the predictor type \n RAtioTimeCutoff: 1e-12; StepsForIncrease: 7;\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	
	std::string::const_iterator iter = configStr.begin();
	std::string::const_iterator end = configStr.end();
	
	
	config::Cauchy<T> structure;
	parsing::ConfigSettingParser<std::string::const_iterator,config::Cauchy<T>, T> parser;
	bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, structure);
	
	
	BOOST_CHECK(parsed && iter == end);
	BOOST_CHECK(fabs(structure.cycle_cutoff_time - 5.76e2) < tol);
	BOOST_CHECK(fabs(structure.ratio_cutoff_time / 1e-12 - 1) < tol);
	
	
}


BOOST_AUTO_TEST_CASE(all_config_settings_d)
{
	using namespace bertini::classic;
	using namespace bertini::tracking;
	double tol = 1e-15;

	SplitInputFile inputfile = parsing::ParseInputFile("Config \n ODEPredictor: 7; \n MPType: 0; \n MaxNewtonIts: 7;  FinalTol: 1.845e-7;\n SampleFactor: 0.647; \n NumSamplePoints: 7;\n\n end;  \n iNpUt % \n  \n variable x; \n ENd;");
	
	
	std::string configStr = inputfile.Config();
	
	auto sets = parsing::GetConfigSettings<double, config::Predictor, config::Newton, config::Tolerances<double>, config::Endgame<double>, config::Stepping<double>>(configStr);
	
	config::Stepping<double> steps = std::get<config::Stepping<double>>(sets);
	config::Predictor pred = std::get<config::Predictor>(sets);
	config::Newton newt = std::get<config::Newton>(sets);
	config::Tolerances<double> tols = std::get<config::Tolerances<double>>(sets);
	config::Endgame<double> end = std::get<config::Endgame<double>>(sets);
	
	BOOST_CHECK(pred == config::Predictor::RKDormandPrince56);
	BOOST_CHECK(fabs(steps.max_step_size - 0.1) < tol);
	BOOST_CHECK(newt.max_num_newton_iterations == 7);
	BOOST_CHECK(newt.min_num_newton_iterations == 1);
	BOOST_CHECK(fabs(tols.final_tolerance - 1.845e-7) < tol);
	BOOST_CHECK(fabs(tols.final_tolerance_multiplier - 10) < tol);
	BOOST_CHECK(fabs(end.sample_factor - 0.647) < tol);
	BOOST_CHECK(end.num_sample_points == 7);
	BOOST_CHECK(fabs(end.min_track_time - 1e-100) < tol);
}





BOOST_AUTO_TEST_SUITE_END()
