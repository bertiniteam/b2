//This file is part of Bertini 2.
//
//bertini2/io/parsing/settings_rules/algorithm.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/setting_rules/algorithm.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/setting_rules/algorithm.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/setting_rules/algorithm.hpp
 
 \brief Provides the parsing rules for algorithm-related settings in bertini2.
 */

#pragma once


#include "bertini2/io/parsing/settings_rules/base.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"

namespace bertini {
	namespace parsing {
		namespace classic {


			
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
				
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "ConfigMidPath")
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