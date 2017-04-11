//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/setting_rules.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/setting_rules.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/setting_rules.hpp
 
 \brief Provides the parsing rules for settings in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/qi_files.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"

#include "bertini2/nag_algorithms/common/config.hpp"

namespace bertini {
	namespace parsing {
		namespace classic {

			using namespace bertini::tracking;
			namespace qi = ::boost::spirit::qi;
			namespace ascii = ::boost::spirit::ascii;
			
			
			/**
			 Qi Parser object for parsing config settings  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
			 This is a base case template, which we specialize below
			 */
			template<typename Iterator, typename Structure, typename T = double, typename Skipper = ascii::space_type> //boost::spirit::unused_type
			struct ConfigSettingParser : qi::grammar<Iterator, Structure(), Skipper>
			{ 
			};
			
			
			
			
			/**
			 Qi Parser object for parsing config settings  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
			 To use this parser, construct an object of its type, then use it to parse.
			 
			 \code
			 ConfigT conf
			 std::string str = "tracktype: 1; tracktolbeforeeg: 1e-8; \n";
			 
			 std::string::const_iterator iter = str.begin();
			 std::string::const_iterator end = str.end();
			 
			 
			 bertini::ConfigSettingParser<std::string::const_iterator, ConfigT> S;
			 
			 
			 bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, conf);
			 
			 \endcode
			 
			 \brief Qi Parser object for parsing config file to determine settings.
			 
			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, PrecisionType, T, Skipper> : qi::grammar<Iterator, PrecisionType(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigPrecisionType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					precisiontype_.add("0", PrecisionType::Fixed);
					precisiontype_.add("1", PrecisionType::Adaptive);
					precisiontype_.add("2", PrecisionType::Adaptive);
					
					
					
					
					root_rule_.name("ConfigPrecisionType_root_rule");
					
					root_rule_ = (config_name_[_val = _1] >> -no_setting_) | no_setting_[_val = PrecisionType::Adaptive];
					
					config_name_.name("precisiontype_");
					config_name_ = *(char_ - (no_case["mptype"]>>':')) >> (no_case["mptype"] >> ':') >> precisiontype_[_val = _1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - no_case["mptype:"]);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, PrecisionType(), ascii::space_type > root_rule_, config_name_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_;
				
				qi::symbols<char,PrecisionType> precisiontype_;
				
				
			}; //re: PrecisionTypeParser
			
			
			
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Predictor, T, Skipper> : qi::grammar<Iterator, config::Predictor(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigPredictorType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					predictor_.add("-1", config::Predictor::Constant);
					predictor_.add("0", config::Predictor::Euler);
					predictor_.add("1", config::Predictor::Heun);
					predictor_.add("2", config::Predictor::RK4);
					predictor_.add("3", config::Predictor::Heun);
					predictor_.add("4", config::Predictor::RKNorsett34);
					predictor_.add("5", config::Predictor::RKF45);
					predictor_.add("6", config::Predictor::RKCashKarp45);
					predictor_.add("7", config::Predictor::RKDormandPrince56);
					predictor_.add("8", config::Predictor::RKVerner67);
					
					
					std::string setting_name = "odepredictor";
					
					
					root_rule_.name("ConfigPredictor_root_rule");
					
					root_rule_ = (config_name_[_val = _1] >> -no_setting_) | no_setting_[_val = config::Predictor::RKF45];
					
					config_name_.name("predictor_");
					config_name_ = *(char_ - (no_case[setting_name] >> ':')) >> (no_case[setting_name] >> ':') >> predictor_[_val = _1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - no_case[setting_name]);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Predictor(), ascii::space_type > root_rule_, config_name_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_;
				
				qi::symbols<char,config::Predictor> predictor_;
				
				
			}; //re: PredictorTypeParser
			
			
			
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Security<T>, T, Skipper> : qi::grammar<Iterator, config::Security<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigSecurityType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string level_name = "securitylevel";
					std::string maxnorm_name = "securitymaxnorm";
					
					
					root_rule_.name("ConfigSecurity_root_rule");
					
					root_rule_ = ((security_level_[phx::bind( [this](config::Security<T> & S, int l)
															 {
																 S.level = l;
															 }, _val, _1 )]
								   ^ security_max_norm_[phx::bind( [this](config::Security<T> & S, T norm)
																  {
																	  S.max_norm = norm;
																  }, _val, _1 )])
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[level_name] >> ':') | (no_case[maxnorm_name] >> ':');
					
					security_level_.name("security_level_");
					security_level_ = *(char_ - all_names_) >> (no_case[level_name] >> ':') >> qi::int_[_val = _1] >> ';';
					
					security_max_norm_.name("security_max_norm_");
					security_max_norm_ = *(char_ - all_names_) >> (no_case[maxnorm_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Security<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, int(), ascii::space_type > security_level_;
				qi::rule<Iterator, T(), ascii::space_type > security_max_norm_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: SecurityParser
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, algorithm::config::Tolerances<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::Tolerances<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigTolerancesType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string newton_before_name = "tracktolbeforeeg";
					std::string newton_during_name = "tracktolduringeg";
					std::string final_tol_name = "finaltol";
					std::string path_trunc_name = "pathtruncationthreshold";
					
					
					root_rule_.name("ConfigTolerances_root_rule");
					
					root_rule_ = ((newton_before_endgame_[phx::bind( [this](algorithm::config::Tolerances<T> & S, T l)
																	{
																		S.newton_before_endgame = l;
																	}, _val, _1 )]
								   ^ newton_during_endgame_[phx::bind( [this](algorithm::config::Tolerances<T> & S, T num)
																	  {
																		  S.newton_during_endgame = num;
																	  }, _val, _1 )]
								   ^ final_tol_[phx::bind( [this](algorithm::config::Tolerances<T> & S, T num)
														  {
															  S.final_tolerance = num;
														  }, _val, _1 )]
								   ^ path_trunc_threshold_[phx::bind( [this](algorithm::config::Tolerances<T> & S, T num)
																	 {
																		 S.path_truncation_threshold = num;
																	 }, _val, _1 )])
								  
								  >> -no_setting_)
					| no_setting_;
					
					
					all_names_ = (no_case[newton_before_name] >> ':') | (no_case[newton_during_name] >> ':')| (no_case[final_tol_name] >> ':')
					| (no_case[path_trunc_name] >> ':');
					
					newton_before_endgame_.name("newton_before_endgame_");
					newton_before_endgame_ = *(char_ - all_names_) >> (no_case[newton_before_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					newton_during_endgame_.name("newton_during_endgame_");
					newton_during_endgame_ = *(char_ - all_names_) >> (no_case[newton_during_name] >>':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					final_tol_.name("final_tol_");
					final_tol_ = *(char_ - all_names_) >> (no_case[final_tol_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					path_trunc_threshold_.name("path_trunc_threshold_");
					path_trunc_threshold_ = *(char_ - all_names_) >> (no_case[path_trunc_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::Tolerances<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > newton_before_endgame_, newton_during_endgame_,
				final_tol_, path_trunc_threshold_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: ToleranceParser
			
			
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Stepping<T>, T, Skipper> : qi::grammar<Iterator, config::Stepping<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigSteppingType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string maxstep_name = "maxstepsize";
					std::string stepsuccess_name = "stepsuccessfactor";
					std::string stepfail_name = "stepfailfactor";
					std::string stepsincrease_name = "stepsforincrease";
					std::string maxnumsteps_name = "maxnumbersteps";
					
					
					root_rule_.name("ConfigStepping_root_rule");
					
					root_rule_ = ((max_step_size_[phx::bind( [this](config::Stepping<T> & S, T num)
															{
																S.max_step_size = num;
															}, _val, _1 )]
								   ^ stepsize_success_[phx::bind( [this](config::Stepping<T> & S, T num)
																 {
																	 S.step_size_success_factor = num;
																 }, _val, _1 )]
								   ^ stepsize_fail_[phx::bind( [this](config::Stepping<T> & S, T num)
															  {
																  S.step_size_fail_factor = num;
															  }, _val, _1 )]
								   ^ steps_increase_[phx::bind( [this](config::Stepping<T> & S, unsigned num)
																	 {
																		 S.consecutive_successful_steps_before_stepsize_increase = num;
																	 }, _val, _1 )]
								   ^ max_num_steps_[phx::bind( [this](config::Stepping<T> & S, unsigned num)
																	{
																		S.max_num_steps = num;
																	}, _val, _1 )])
								  
								  >> -no_setting_)
					| no_setting_;
					
					
					all_names_ = (no_case[maxstep_name] >> ':') | (no_case[stepsuccess_name] >> ':')|
					(no_case[stepfail_name] >> ':')| (no_case[stepsincrease_name] >> ':') | (no_case[maxnumsteps_name] >> ':');
					
					max_step_size_.name("max_step_size_");
					max_step_size_ = *(char_ - all_names_) >> (no_case[maxstep_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					stepsize_success_.name("stepsize_success_");
					stepsize_success_ = *(char_ - all_names_) >> (no_case[stepsuccess_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					stepsize_fail_.name("stepsize_fail_");
					stepsize_fail_ = *(char_ - all_names_) >> (no_case[stepfail_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					steps_increase_.name("steps_increase_");
					steps_increase_ = *(char_ - all_names_) >> (no_case[stepsincrease_name] >> ':') >> qi::uint_[_val=_1] >> ';';
					
					max_num_steps_.name("max_num_steps_");
					max_num_steps_ = *(char_ - all_names_) >> (no_case[maxnumsteps_name] >> ':') >> qi::uint_[_val=_1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Stepping<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > max_step_size_, stepsize_success_,
				stepsize_fail_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > steps_increase_, max_num_steps_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: SteppingParser
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Newton, T, Skipper> : qi::grammar<Iterator, config::Newton(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigNewtonType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string maxits_name = "maxnewtonits";
					
					
					root_rule_.name("ConfigNewton_root_rule");
					
					root_rule_ = (max_its_[phx::bind( [this](config::Newton & S, unsigned num)
													 {
														 S.max_num_newton_iterations = num;
													 }, _val, _1 )]
								  
								  >> -no_setting_)
					| no_setting_;
					
					
					
					all_names_ = eps >> (no_case[maxits_name] >> ':');
					
					max_its_.name("max_its_");
					max_its_ = *(char_ - no_case[maxits_name]) >> (no_case[maxits_name] >> ':') >> qi::uint_[_val=_1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Newton(), Skipper> root_rule_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > max_its_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
			}; //re: NewtonParser
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Endgame<T>, T, Skipper> : qi::grammar<Iterator, config::Endgame<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigEndgameType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string samplefactor_name = "samplefactor";
					std::string numpoints_name = "numsamplepoints";
					std::string mintrack_name = "nbhdradius";
					
					
					root_rule_.name("ConfigEndgame_root_rule");
					
					root_rule_ = ((sample_factor_[phx::bind( [this](config::Endgame<T> & S, T num)
															{
																S.sample_factor = num;
															}, _val, _1 )]
								   ^ min_track_[phx::bind( [this](config::Endgame<T> & S, T num)
														  {
															  S.min_track_time = num;
														  }, _val, _1 )]
								   ^ num_sample_[phx::bind( [this](config::Endgame<T> & S, unsigned num)
															  {
																  S.num_sample_points = num;
															  }, _val, _1 )])
								  
								  >> -no_setting_)
					| no_setting_;
					
					
					all_names_ = (no_case[samplefactor_name] >> ':') | (no_case[numpoints_name] >> ':')| (no_case[mintrack_name] >> ':');
					
					sample_factor_.name("sample_factor_");
					sample_factor_ = *(char_ - all_names_) >> (no_case[samplefactor_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					min_track_.name("min_track_");
					min_track_ = *(char_ - all_names_) >> (no_case[mintrack_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					num_sample_.name("num_sample_");
					num_sample_ = *(char_ - all_names_) >> (no_case[numpoints_name] >> ':')
					>> qi::uint_[_val=_1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Endgame<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > sample_factor_, min_track_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > num_sample_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
			}; //re: EndgameParser
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::PowerSeries, T, Skipper> : qi::grammar<Iterator, config::PowerSeries(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigPowerSeriesType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string maxcycle_name = "maxcyclenum";
					
					
					root_rule_.name("ConfigPowerSeries_root_rule");
					
					root_rule_ = (max_cycle_[phx::bind( [this](config::PowerSeries & S, unsigned num)
													   {
														   S.max_cycle_number = num;
													   }, _val, _1 )]
								  
								  >> -no_setting_)
					| no_setting_;
					
					
					
					all_names_ = eps >> (no_case[maxcycle_name] >> ':');
					
					max_cycle_.name("max_cycle_");
					max_cycle_ = *(char_ - all_names_) >> (no_case[maxcycle_name] >> ':') >> qi::uint_[_val=_1] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::PowerSeries(), Skipper> root_rule_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > max_cycle_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
			}; //re: PowerSeriesParser
			
			
			
			
			
			/**

			 */
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::Cauchy<T>, T, Skipper> : qi::grammar<Iterator, config::Cauchy<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigCauchyType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string cyclecutoff_name = "cycletimecutoff";
					std::string ratiocutoff_name = "ratiotimecutoff";
					
					
					root_rule_.name("ConfigCauchy_root_rule");
					
					root_rule_ = ((cycle_cutoff_[phx::bind( [this](config::Cauchy<T> & S, T num)
														   {
															   S.cycle_cutoff_time = num;
														   }, _val, _1 )]
								   ^ ratio_cutoff_[phx::bind( [this](config::Cauchy<T> & S, T num)
															 {
																 S.ratio_cutoff_time = num;
															 }, _val, _1 )])
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[cyclecutoff_name] >> ':') | (no_case[ratiocutoff_name] >> ':');
					
					cycle_cutoff_.name("cycle_cutoff_");
					cycle_cutoff_ = *(char_ - all_names_) >> (no_case[cyclecutoff_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					ratio_cutoff_.name("ratio_cutoff_");
					ratio_cutoff_ = *(char_ - all_names_) >> (no_case[ratiocutoff_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( [this](T & num, std::string str)
														   {
															   num = bertini::NumTraits<T>::FromString(str);
														   }, _val, _1 )] >> ';';
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::Cauchy<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > cycle_cutoff_, ratio_cutoff_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: CauchyParser

			

			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, config::TrackBack, T, Skipper> : qi::grammar<Iterator, config::TrackBack(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigTrackBack")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string min_cycle_name = "mincycletrackback";
					std::string junk_removal_name = "junkremovaltrackback";
					std::string maxdepth_name = "maxldtdepth";

					root_rule_.name("ConfigTrackBack_root_rule");
					
					root_rule_ = root_rule_ = ((min_cycle_[phx::bind( [this](config::TrackBack & S, unsigned num)
																	   {
																		   S.minimum_cycle = num;
																	   }, _val, _1 )]
											   ^ junk_removal_[phx::bind( [this](config::TrackBack & S, int num)
																		 {
																			 S.junk_removal_test = static_cast<bool>(num);
																		 }, _val, _1 )]
											   ^ max_depth_[phx::bind( [this](config::TrackBack & S, unsigned num)
																	   {
																		   S.max_depth_LDT = num;
																	   }, _val, _1 )]
											  )
											  >> -no_setting_) | no_setting_;
					
					
					
					all_names_ = eps >> (no_case[min_cycle_name] >> ':') |
										(no_case[junk_removal_name] >> ':') |
										(no_case[maxdepth_name] >> ':');
					
					min_cycle_.name("min_cycle_");
					min_cycle_ = *(char_ - all_names_) >> (no_case[min_cycle_name] >> ':') >> qi::uint_[_val=_1] >> ';';
					
					junk_removal_.name("junk_removal_");
					junk_removal_ = *(char_ - all_names_) >> (no_case[junk_removal_name] >> ':') >> qi::int_[_val=_1] >> ';';

					max_depth_.name("max_depth_");
					max_depth_ = *(char_ - all_names_) >> (no_case[maxdepth_name] >> ':') >> qi::uint_[_val=_1] >> ';';

					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, config::TrackBack(), Skipper> root_rule_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > min_cycle_, max_depth_;
				qi::rule<Iterator, int(), ascii::space_type > junk_removal_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
			}; //re: TrackBackParser


			/**
			As of this writing, there are no meaningful settings in FixedPrecisionConfig, 
			so this parser is ... empty
			*/
			template<typename Iterator, typename T, typename Skipper>
			struct ConfigSettingParser<Iterator, config::FixedPrecisionConfig<T>, T, Skipper> : qi::grammar<Iterator, config::FixedPrecisionConfig<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigFixedPrecision")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					

					
					root_rule_.name("ConfigFixedPrecision_root_rule");
					root_rule_ = eps[_val = config::FixedPrecisionConfig<T>()] >> omit[*(char_)];
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
				}
				
				
			private:
				qi::rule<Iterator, config::FixedPrecisionConfig<T>(), ascii::space_type > root_rule_;
			}; //re: FixedPrecisionConfig Parser

			





			using AdaptiveMultiplePrecisionConfig = config::AdaptiveMultiplePrecisionConfig;
			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, AdaptiveMultiplePrecisionConfig, T, Skipper> : qi::grammar<Iterator, AdaptiveMultiplePrecisionConfig(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "AdaptiveMultiplePrecisionConfig")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string coefficient_bound_name = "coefficientbound";
					std::string degree_bound_name = "degreebound";
					std::string lin_solve_error_bnd_name = "epsilon";
					std::string jac_eval_err_bnd_name = "phi";
					std::string func_eval_err_bnd_name = "psi";
					std::string safety_one_name = "ampsafetydigits1";
					std::string safety_two_name = "ampsafetydigits2";
					std::string max_prec_name = "ampmaxprec";
					std::string consec_steps_prec_dec_name = "maxstepsprecisiondecrease";
					std::string max_num_prec_decs_name = "maxnumprecdecreases";

					root_rule_.name("ConfigAMP_root_rule");
					
					root_rule_ = ((coefficient_bound_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, mpfr_float num)
														   {
															   S.coefficient_bound = num;
														   }, _val, _1 )]
								   ^ degree_bound_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, mpfr_float num)
															 {
																 S.degree_bound = num;
															 }, _val, _1 )]
								   ^ lin_solve_error_bnd_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, mpfr_float num)
															 {
																 S.epsilon = num;
															 }, _val, _1 )]
								   ^ jac_eval_err_bnd_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, mpfr_float num)
															 {
																 S.Phi = num;
															 }, _val, _1 )]
								   ^ func_eval_err_bnd_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, mpfr_float num)
															 {
																 S.Psi = num;
															 }, _val, _1 )]
								   ^ safety_one_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, int num)
															 {
																 S.safety_digits_1 = num;
															 }, _val, _1 )]
								   ^ safety_two_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, int num)
															 {
																 S.safety_digits_2 = num;
															 }, _val, _1 )]
								   ^ max_prec_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, unsigned num)
															 {
																 S.maximum_precision = num;
															 }, _val, _1 )]
								   ^ consec_steps_prec_dec_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, unsigned num)
															 {
																 S.consecutive_successful_steps_before_precision_decrease = num;
															 }, _val, _1 )]
								   ^ max_num_prec_decs_[phx::bind( [this](config::AdaptiveMultiplePrecisionConfig & S, unsigned num)
															 {
																 S.max_num_precision_decreases = num;
															 }, _val, _1 )]
								   )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[coefficient_bound_name] >> ':') | 
								 (no_case[degree_bound_name] >> ':') |
								 (no_case[lin_solve_error_bnd_name] >> ':') |
								 (no_case[jac_eval_err_bnd_name] >> ':') |
								 (no_case[func_eval_err_bnd_name] >> ':') |
								 (no_case[safety_one_name] >> ':') |
								 (no_case[safety_two_name] >> ':') |
								 (no_case[max_prec_name] >> ':') |
								 (no_case[consec_steps_prec_dec_name] >> ':') |
								 (no_case[max_num_prec_decs_name] >> ':')
								 ;
					

					auto str_to_T = [this](mpfr_float & num, std::string str)
									   {
									   	std::cout << str << std::endl;
										   num = bertini::NumTraits<T>::FromString(str);
									   };


					degree_bound_.name("degree_bound_");
					degree_bound_ = *(char_ - all_names_) >> (no_case[degree_bound_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';
					
					coefficient_bound_.name("coefficient_bound_");
					coefficient_bound_ = *(char_ - all_names_) >> (no_case[coefficient_bound_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';
					
					lin_solve_error_bnd_.name("lin_solve_error_bnd_");
					lin_solve_error_bnd_ = *(char_ - all_names_) >> (no_case[lin_solve_error_bnd_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					jac_eval_err_bnd_.name("jac_eval_err_bnd_");
					jac_eval_err_bnd_ = *(char_ - all_names_) >> (no_case[jac_eval_err_bnd_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					func_eval_err_bnd_.name("func_eval_err_bnd_");
					func_eval_err_bnd_ = *(char_ - all_names_) >> (no_case[func_eval_err_bnd_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					safety_one_.name("safety_one_");
					safety_one_ = *(char_ - all_names_) >> (no_case[safety_one_name] >> ':') >> qi::int_[_val=_1] >> ';';

					safety_two_.name("safety_two_");
					safety_two_ = *(char_ - all_names_) >> (no_case[safety_two_name] >> ':') >> qi::int_[_val=_1] >> ';';

					max_prec_.name("max_prec_");
					max_prec_ = *(char_ - all_names_) >> (no_case[max_prec_name] >> ':') >> qi::uint_[_val=_1] >> ';';

					consec_steps_prec_dec_.name("consec_steps_prec_dec_");
					consec_steps_prec_dec_ = *(char_ - all_names_) >> (no_case[consec_steps_prec_dec_name] >> ':') >> qi::uint_[_val=_1] >> ';';

					max_num_prec_decs_.name("max_num_prec_decs_");
					max_num_prec_decs_ = *(char_ - all_names_) >> (no_case[max_num_prec_decs_name] >> ':') >> qi::uint_[_val=_1] >> ';';


					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, AdaptiveMultiplePrecisionConfig(), ascii::space_type > root_rule_;

				qi::rule<Iterator, mpfr_float(), ascii::space_type > degree_bound_, coefficient_bound_, lin_solve_error_bnd_, jac_eval_err_bnd_, func_eval_err_bnd_;
				qi::rule<Iterator, int(), ascii::space_type > safety_one_, safety_two_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > max_prec_, consec_steps_prec_dec_, max_num_prec_decs_;

				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: AMPParser



			template<typename Iterator, typename T, typename Skipper> 
			struct ConfigSettingParser<Iterator, algorithm::config::ZeroDim<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::ZeroDim<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "algorithm::config::ZeroDim")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string init_prec_name = "initialprec";
					std::string max_cross_resolve_name = "maxcrossedpathresolves";
					std::string start_time_name = "starttime";
					std::string endgame_boundary_name = "endgamebdry";
					std::string target_time_name = "targettime";
					std::string path_variable_name_name = "pathvarname";


					root_rule_.name("ConfigZeroDim_root_rule");
					
					root_rule_ = ((init_prec_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, int num)
														   {
															   S.initial_ambient_precision = num;
														   }, _val, _1 )]
								   ^ path_variable_name_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, std::string omnom)
															 {
																 S.path_variable_name = omnom;
															 }, _val, _1 )]
								   ^ max_cross_resolve_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, int num)
															 {
																 S.max_num_crossed_path_resolve_attempts = num;
															 }, _val, _1 )]
								   ^ start_time_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, T num)
															 {
																 S.start_time = num;
															 }, _val, _1 )]
								   ^ endgame_boundary_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, T num)
															 {
																 S.endgame_boundary = num;
															 }, _val, _1 )]
								   ^ target_time_[phx::bind( [this](algorithm::config::ZeroDim<T> & S, T num)
															 {
																 S.target_time = num;
															 }, _val, _1 )]
								   )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[init_prec_name] >> ':') | 
								 (no_case[target_time_name] >> ':') |
								 (no_case[path_variable_name_name] >> ':') |
								 (no_case[max_cross_resolve_name] >> ':') |
								 (no_case[start_time_name] >> ':') |
								 (no_case[endgame_boundary_name] >> ':')
								 ;
					

					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };



					valid_variable_name_.name("valid_variable_name_");
					valid_variable_name_ = +qi::alpha >> *(qi::alnum | qi::char_("[]_") );


					path_variable_name_.name("path_variable_name_");
					path_variable_name_ = *(char_ - all_names_) >> (no_case[path_variable_name_name] >> ':') >> valid_variable_name_ >> ';';

					init_prec_.name("init_prec_");
					init_prec_ = *(char_ - all_names_) >> (no_case[init_prec_name] >> ':') >> qi::int_[_val=_1] >> ';';

					max_cross_resolve_.name("max_cross_resolve_");
					max_cross_resolve_ = *(char_ - all_names_) >> (no_case[max_cross_resolve_name] >> ':') >> qi::int_[_val=_1] >> ';';

					start_time_.name("start_time_");
					start_time_ = *(char_ - all_names_) >> (no_case[start_time_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					endgame_boundary_.name("endgame_boundary_");
					endgame_boundary_ = *(char_ - all_names_) >> (no_case[endgame_boundary_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					target_time_.name("target_time_");
					target_time_ = *(char_ - all_names_) >> (no_case[target_time_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';



					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::ZeroDim<T>(), ascii::space_type > root_rule_;

				qi::rule<Iterator, T(), ascii::space_type > start_time_, target_time_, endgame_boundary_;
				qi::rule<Iterator, int(), ascii::space_type > init_prec_, max_cross_resolve_;
				qi::rule<Iterator, std::string(), ascii::space_type > valid_variable_name_, path_variable_name_;
				

				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: ZeroDim

			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, algorithm::config::MidPath<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::MidPath<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigCauchyType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string same_point_tol_name = "midpathtol";
					
					
					root_rule_.name("ConfigMidPath_root_rule");
					
					root_rule_ = ((same_point_tol_[phx::bind( [this](algorithm::config::MidPath<T> & S, T num)
														   {
															   S.same_point_tolerance = num;
														   }, _val, _1 )]
								  )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[same_point_tol_name] >> ':');
					
					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };


					same_point_tol_.name("same_point_tol_");
					same_point_tol_ = *(char_ - all_names_) >> (no_case[same_point_tol_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';
					
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::MidPath<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > same_point_tol_;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
			}; //re: MidPathParser


			template<typename Iterator, typename T, typename Skipper> //boost::spirit::unused_type
			struct ConfigSettingParser<Iterator, algorithm::config::AutoRetrack<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::AutoRetrack<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigAutoretrack")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string decrease_factor_name = "retracktolfactor";
					
					
					root_rule_.name("ConfigAutoRetrack_root_rule");
					
					root_rule_ = ((decrease_factor[phx::bind( [this](algorithm::config::AutoRetrack<T> & S, T num)
														   {
															   S.midpath_decrease_tolerance_factor = num;
														   }, _val, _1 )]
								  )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[decrease_factor_name] >> ':');
					
					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };


					decrease_factor.name("decrease_factor");
					decrease_factor = *(char_ - all_names_) >> (no_case[decrease_factor_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';
					
					
					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::AutoRetrack<T>(), ascii::space_type > root_rule_;
				qi::rule<Iterator, T(), ascii::space_type > decrease_factor;
				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
			}; //re: AutoRetrackParser


			template<typename Iterator, typename T, typename Skipper> 
			struct ConfigSettingParser<Iterator, algorithm::config::Sharpening<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::Sharpening<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "algorithm::config::Sharpening")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string sharpen_digits_name = "sharpendigits";
					std::string func_res_tol_name = "functiontolerance";
					std::string ratio_tol_name = "ratiotolerance";


					root_rule_.name("ConfigSharpening_root_rule");
					
					root_rule_ = ((sharpen_digits_[phx::bind( [this](algorithm::config::Sharpening<T> & S, unsigned num)
														   {
															   S.sharpendigits = num;
														   }, _val, _1 )]
								   ^ func_res_tol_[phx::bind( [this](algorithm::config::Sharpening<T> & S, T num)
															 {
																 S.function_residual_tolerance = num;
															 }, _val, _1 )]
								   ^ ratio_tol_[phx::bind( [this](algorithm::config::Sharpening<T> & S, T num)
															 {
																 S.ratio_tolerance = num;
															 }, _val, _1 )]
								   )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[sharpen_digits_name] >> ':') | 
								 (no_case[func_res_tol_name] >> ':') |
								 (no_case[ratio_tol_name] >> ':')
								 ;
					

					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };


					func_res_tol_.name("func_res_tol_");
					func_res_tol_ = *(char_ - all_names_) >> (no_case[func_res_tol_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';
					
					ratio_tol_.name("ratio_tol_");
					ratio_tol_ = *(char_ - all_names_) >> (no_case[ratio_tol_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					sharpen_digits_.name("sharpen_digits_");
					sharpen_digits_ = *(char_ - all_names_) >> (no_case[sharpen_digits_name] >> ':') >> qi::uint_[_val=_1] >> ';';



					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::Sharpening<T>(), ascii::space_type > root_rule_;

				qi::rule<Iterator, T(), ascii::space_type > func_res_tol_, ratio_tol_;
				qi::rule<Iterator, unsigned int(), ascii::space_type > sharpen_digits_;

				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: Sharpening

			template<typename Iterator, typename T, typename Skipper> 
			struct ConfigSettingParser<Iterator, algorithm::config::Regeneration<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::Regeneration<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "algorithm::config::Regeneration")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string regen_remove_inf_name = "regenremoveinf";
					std::string slice_before_name = "slicetolbeforeeg";
					std::string slice_during_name = "slicetolduringeg";
					std::string slice_final_name = "slicefinaltol";
					std::string higher_dim_check_name = "regenhigherdimtest";
					std::string start_level_name = "regenstartlevel";


					root_rule_.name("ConfigRegeneration_root_rule");
					
					root_rule_ = ((regen_remove_inf_[phx::bind( [this](algorithm::config::Regeneration<T> & S, int num)
														   {
															   S.remove_infinite_endpoints = static_cast<bool>(num);
														   }, _val, _1 )]
								   ^ higher_dim_check_[phx::bind( [this](algorithm::config::Regeneration<T> & S, int num)
															 {
																 S.higher_dimension_check = static_cast<bool>(num);
															 }, _val, _1 )]
								   ^ start_level_[phx::bind( [this](algorithm::config::Regeneration<T> & S, int num)
															 {
																 S.start_level = num;
															 }, _val, _1 )]
								   ^ slice_before_[phx::bind( [this](algorithm::config::Regeneration<T> & S, T num)
															 {
																 S.newton_before_endgame = num;
															 }, _val, _1 )]
								   ^ slice_during_[phx::bind( [this](algorithm::config::Regeneration<T> & S, T num)
															 {
																 S.newton_during_endgame = num;
															 }, _val, _1 )]
								   ^ slice_final_[phx::bind( [this](algorithm::config::Regeneration<T> & S, T num)
															 {
																 S.final_tolerance = num;
															 }, _val, _1 )]
								   )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[regen_remove_inf_name] >> ':') | 
								 (no_case[higher_dim_check_name] >> ':') |
								 (no_case[start_level_name] >> ':') |
								 (no_case[slice_before_name] >> ':') |
								 (no_case[slice_during_name] >> ':') |
								 (no_case[slice_final_name] >> ':')
								 ;
					

					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };


					higher_dim_check_.name("higher_dim_check_");
					higher_dim_check_ = *(char_ - all_names_) >> (no_case[higher_dim_check_name] >> ':') >> qi::int_[_val=_1] >> ';';

					start_level_.name("start_level_");
					start_level_ = *(char_ - all_names_) >> (no_case[start_level_name] >> ':') >> qi::int_[_val=_1] >> ';';

					regen_remove_inf_.name("regen_remove_inf_");
					regen_remove_inf_ = *(char_ - all_names_) >> (no_case[regen_remove_inf_name] >> ':') >> qi::int_[_val=_1] >> ';';



					slice_before_.name("slice_before_");
					slice_before_ = *(char_ - all_names_) >> (no_case[slice_before_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					slice_during_.name("slice_during_");
					slice_during_ = *(char_ - all_names_) >> (no_case[slice_during_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					slice_final_.name("slice_final_");
					slice_final_ = *(char_ - all_names_) >> (no_case[slice_final_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';



					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::Regeneration<T>(), ascii::space_type > root_rule_;

				qi::rule<Iterator, T(), ascii::space_type > slice_before_, slice_during_, slice_final_;
				qi::rule<Iterator, int(), ascii::space_type > regen_remove_inf_, start_level_, higher_dim_check_;

				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: Regeneration


			template<typename Iterator, typename T, typename Skipper> 
			struct ConfigSettingParser<Iterator, algorithm::config::PostProcessing<T>, T, Skipper> : qi::grammar<Iterator, algorithm::config::PostProcessing<T>(), Skipper>
			{
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "algorithm::config::PostProcessing")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					
					
					std::string real_thresh_name = "imagthreshold";
					std::string endpoint_finite_name = "endpointfinitethreshold";
					std::string same_point_name = "endpointsamethreshold";



					root_rule_.name("ConfigPostProcessing_root_rule");
					
					root_rule_ = ((real_threshold_[phx::bind( [this](algorithm::config::PostProcessing<T> & S, T num)
															 {
																 S.real_threshold = num;
															 }, _val, _1 )]
								   ^ endpoint_finite_[phx::bind( [this](algorithm::config::PostProcessing<T> & S, T num)
															 {
																 S.endpoint_finite_threshold = num;
															 }, _val, _1 )]
								   ^ same_point_[phx::bind( [this](algorithm::config::PostProcessing<T> & S, T num)
															 {
																 S.same_point_tolerance = num;
															 }, _val, _1 )]
								   )
								  >> -no_setting_) | no_setting_;
					
					all_names_ = (no_case[real_thresh_name] >> ':') |
								 (no_case[endpoint_finite_name] >> ':') |
								 (no_case[same_point_name] >> ':')
								 ;
					

					auto str_to_T = [this](T & num, std::string str)
									   {
										   num = bertini::NumTraits<T>::FromString(str);
									   };

					real_threshold_.name("real_thresh_");
					real_threshold_ = *(char_ - all_names_) >> (no_case[real_thresh_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					endpoint_finite_.name("endpoint_finite_");
					endpoint_finite_ = *(char_ - all_names_) >> (no_case[endpoint_finite_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';

					same_point_.name("same_point_");
					same_point_ = *(char_ - all_names_) >> (no_case[same_point_name] >> ':')
					>> mpfr_rules.number_string_[phx::bind( str_to_T, _val, _1 )] >> ';';



					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
					
					no_decl_.name("no_decl_");
					no_decl_ = *(char_);
					
					
					
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
					
				}
				
				
			private:
				qi::rule<Iterator, algorithm::config::PostProcessing<T>(), ascii::space_type > root_rule_;

				qi::rule<Iterator, T(), ascii::space_type > real_threshold_, endpoint_finite_, same_point_;

				qi::rule<Iterator, ascii::space_type, std::string()> no_decl_, no_setting_, all_names_;
				MPParserRules<Iterator> mpfr_rules;
				
				
				
			}; //re: PostProcessing




		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini