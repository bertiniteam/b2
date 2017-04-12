//This file is part of Bertini 2.
//
//bertini2/io/parsing/settings_parsers/tracking.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/settings_parsers/tracking.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/settings_parsers/tracking.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/settings_parsers/tracking.hpp
 
 \brief Provides the parsing rules for tracking-related settings in bertini2.
 */

#pragma once


#include "bertini2/io/parsing/settings_parsers/base.hpp"
#include "bertini2/trackers/config.hpp"


namespace bertini {
	namespace parsing {
		namespace classic {

			
			namespace {
				using namespace bertini::tracking;
			}
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
				rules::LongNum<Iterator> mpfr_rules;
				
				
				
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
			template<typename Iterator, typename T, typename Skipper> 
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
				rules::LongNum<Iterator> mpfr_rules;
				
				
				
			}; //re: AMPParser





		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini