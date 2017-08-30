//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/system_rules.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/system_rules.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/system_rules.hpp
 
 \brief Provides the parsing rules for systems in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/qi_files.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/function_rules.hpp"



namespace bertini {
	namespace parsing {
		namespace classic {
			// a few local using statements to reduce typing etc.
			using Function = node::Function;
			using Variable = node::Variable;
			using Node = node::Node;
			
			
			using Fn = std::shared_ptr<Function>;
			using Var = std::shared_ptr<Variable>;
			using Nd = std::shared_ptr<Node>;
			
			/**
			 Qi Parser object for parsing text into the System class.  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
			 To use this parser, construct an object of its type, then use it to parse.
			 
			 \code
			 System sys;
			 std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";
			 
			 std::string::const_iterator iter = str.begin();
			 std::string::const_iterator end = str.end();
			 
			 
			 bertini::SystemParser<std::string::const_iterator> S;
			 
			 
			 bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);
			 
			 \endcode
			 
			 \brief Qi Parser object for parsing text into the System class.
			 
			 This parser could not have been written without the help of SO user sehe.
			 */
			template<typename Iterator, typename Skipper = ascii::space_type>
			struct SystemParser : qi::grammar<Iterator, System(), Skipper>
			{
				
				
				SystemParser() :	function_parser_(&encountered_symbols_),
				/*initialize here with address of encountered_symbols*/
				SystemParser::base_type(root_rule_)
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using boost::spirit::lexeme;
					
					
					declarative_symbols_.add("variable_group",0);
					declarative_symbols_.add("hom_variable_group",1);
					declarative_symbols_.add("variable",2);
					declarative_symbols_.add("function",3);
					declarative_symbols_.add("constant",4);
					declarative_symbols_.add("parameter",5);
					declarative_symbols_.add("implicit_parameter",6);
					declarative_symbols_.add("pathvariable",7);
					declarative_symbols_.add("random",8);
					declarative_symbols_.add("random_real",9);
					
					
					
					
					
					
					special_numbers_.add("pi", node::Pi());
					special_numbers_.add("Pi", special_numbers_.at("pi"));
					
					special_numbers_.add("e", node::E());
					special_numbers_.add("E", special_numbers_.at("e"));
					
					
					special_numbers_.add("i", node::I());
					special_numbers_.add("I", special_numbers_.at("i"));
					
					
					
					encountered_symbols_.add("pi", special_numbers_.at("pi"));
					encountered_symbols_.add("Pi", special_numbers_.at("pi"));
					
					encountered_symbols_.add("e", special_numbers_.at("e"));
					encountered_symbols_.add("E", special_numbers_.at("e"));
					
					
					encountered_symbols_.add("i", special_numbers_.at("i"));
					encountered_symbols_.add("I", special_numbers_.at("i"));
					
					
					
					
					
					//TODO refine this so that counts are enforced at parse time?
					root_rule_.name("system_parsing");
					root_rule_ =
					*(
					  variable_group_ [phx::bind(&System::AddVariableGroup, _val, _1)]
					  |
					  hom_variable_group_ [phx::bind(&System::AddHomVariableGroup, _val, _1)]
					  |
					  variables_ [phx::bind(&System::AddUngroupedVariables, _val, _1)]
					  |
					  functions_ [phx::bind(&System::AddFunctions, _val, _1)]
					  |
					  constants_ [phx::bind(&System::AddConstants, _val, _1)]
					  |
					  parameters_ [phx::bind(&System::AddParameters, _val, _1)]
					  |
					  implicit_parameters_ [phx::bind(&System::AddImplicitParameters, _val, _1)]
					  |
					  path_variable_ [phx::bind(&System::AddPathVariable, _val, _1)]
					  |
					  subfunction_ [phx::bind(&System::AddSubfunction, _val, _1)]
					  |
					  definition_
					  )
					;
					
					
					
					
					
					variables_.name("variables_"); hom_variable_group_.name("hom_variable_group_"); variable_group_.name("variable_group_"); implicit_parameters_.name("implicit_parameters_");
					
					
					variables_			= "variable" > genericvargp_ > ';';
					hom_variable_group_ = "hom_variable_group" > genericvargp_ > ';';
					variable_group_		= "variable_group" > genericvargp_ > ';';
					implicit_parameters_ = "implicit_parameter" > genericvargp_ > ';';
					
					
					
					path_variable_.name("path_variable_");
					path_variable_ = "pathvariable" > new_variable_ > ';';
					
					
					
					genericvargp_.name("genericvargp_");
					genericvargp_ = new_variable_ % ',';
					
					new_variable_.name("new_variable_");
					new_variable_ = unencountered_symbol_ [boost::phoenix::bind( [this](Var & V, std::string str)
																				{
																					MakeAndAddVariable(V,str);
																				}, _val, _1 )];
					
					
					
					functions_.name("functions_"); constants_.name("constants_"); parameters_.name("parameters_");
					
					functions_ = "function" > genericfuncgp_ > ';';
					constants_ = "constant" > genericfuncgp_ > ';';
					parameters_ = "parameter" > genericfuncgp_ > ';';
					
					
					genericfuncgp_.name("genericfuncgp_");
					genericfuncgp_ = new_function_ % ',';
					
					
					new_function_.name("new_function_");
					new_function_ = unencountered_symbol_ [boost::phoenix::bind( [this](Fn & F, std::string str)
																				{
																					MakeAndAddFunction(F,str);
																				}, _val, _1 )];//[_val = make_shared_<Function>() (_1)]
					
					
					
					
					// this rule gets a string.
					unencountered_symbol_.name("unencountered_symbol_");
					unencountered_symbol_ = valid_variable_name_ - lexeme[( declarative_symbols_ | encountered_symbols_ )];
					// i am unsure about the use of lexeme in the above rule (unencountered_symbol).
					
					
					
					
					// get a string which fits the naming rules.
					valid_variable_name_.name("valid_variable_name_");
					valid_variable_name_ = +qi::alpha >> *(qi::alnum | qi::char_("[]_") );
					
					
					
					
					definition_.name("definition_");
					definition_ = (encountered_functions_ > '=' > function_parser_ > ';') [phx::bind( [](const Fn & F, const Nd & N)
																									 {
																										 F->SetRoot(N);
																									 },_1, _2)] ;
					
					
					using qi::_a;
					using qi::omit;
					subfunction_.name("subfunction");
					subfunction_ = new_function_ [_a = _1]  > '=' >
					function_parser_ [_val = _a, phx::bind( [](Fn & F, const Nd & N)
														   {
															   F->SetRoot(N);
														   },_a, _1)]
					// omit close
					> ';';
					
					
					//			debug(root_rule_);
					//
					//
					//
					//			debug(functions_);
					//			debug(constants_);
					//			debug(parameters_);
					//
					//			debug(genericfuncgp_);
					//			debug(new_function_);
					//
					//
					//
					//			debug(definition_);
					//
					//			debug(subfunction_);
					//
					//
					//
					//
					//			debug(variables_);
					//			debug(hom_variable_group_);
					//			debug(variable_group_);
					//			debug(implicit_parameters_);
					//			debug(path_variable_);
					//
					//			debug(new_variable_);
					//			debug(genericvargp_); debug(variable_group_);
					//
					//			debug(unencountered_symbol_);
					//
					//			debug(valid_variable_name_);
					
					
					
					
					//			BOOST_SPIRIT_DEBUG_NODES( (unencountered_symbol_) (new_variable_) (genericvargp_))
					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("System parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
					
				}
				
				
				
				
			private:
				
				// rule declarations.  these are member variables for the parser.
				qi::rule<Iterator, System(), Skipper > root_rule_;
				
				
				qi::rule<Iterator, VariableGroup(), Skipper > variable_group_, hom_variable_group_, variables_, implicit_parameters_;
				qi::rule<Iterator, VariableGroup(), Skipper > genericvargp_;
				
				qi::rule<Iterator, Var(), Skipper> path_variable_;
				qi::rule<Iterator, Var()> new_variable_;
				
				
				
				
				qi::rule<Iterator, std::vector<Fn>(), Skipper > functions_, constants_, parameters_;
				qi::rule<Iterator, std::vector<Fn>(), Skipper > genericfuncgp_;
				qi::rule<Iterator, Fn(), Skipper, qi::locals<Fn> >  subfunction_;
				
				qi::rule<Iterator, Fn()>  new_function_;
				
				
				qi::rule<Iterator, std::string()> unencountered_symbol_;
				
				
				// the rule which determines valid variable names
				qi::rule<Iterator, std::string()> valid_variable_name_;
				
				qi::rule<Iterator, Skipper, qi::unused_type> definition_;
				
				// symbol declarations
				qi::symbols<char,Nd> encountered_symbols_;
				qi::symbols<char,int> declarative_symbols_;
				qi::symbols<char,Fn>  encountered_functions_;
				qi::symbols<char,Nd> special_numbers_;
				
				FunctionParser<Iterator> function_parser_;
				
				/**
				 To accompany the rule for making new functions when you encounter a new symbol.
				 Simultaneously makes a new function, and adds it to the set of symbols.
				 */
				void MakeAndAddFunction(Fn & F, std::string str)
				{
					F = std::make_shared<Function>(str);
					encountered_symbols_.add(str, F);
					encountered_functions_.add(str,F);
				}
				
				/**
				 To accompany the rule for making new variables when you encounter a new symbol.
				 Simultaneously makes a new variable, and adds it to the set of symbols.
				 */
				void MakeAndAddVariable(Var & V, std::string str)
				{
					V = MakeVariable(str);
					encountered_symbols_.add(str, V);
				}
				
				
				void SetRootNode(Fn & F, const Nd & N)
				{
					F->SetRoot(N);
				}
			};
			
		} // re: namespace classic

	} //re: namespace parsing
	
	
	

} // re: namespace bertini
