//This file is part of Bertini 2.0.
//
//system_parsing.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system_parsing.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system_parsing.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Daniel Brake on 2015 June 08.
//
//
// system_parsing.hpp:  This file contains the parser to parse systems from Bertini Classic input files.



#ifndef BERTINI_SYSTEM_PARSING_HPP
#define BERTINI_SYSTEM_PARSING_HPP



#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>




#include "system.hpp"
#include "function_tree/function_parsing.hpp"




#define BOOST_SPIRIT_USE_PHOENIX_V3 1


#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <iostream>
#include <string>















namespace bertini {
	
	// a few local using statements to reduce typing etc.
	using Fn = std::shared_ptr<Function>;
	using Var = std::shared_ptr<Variable>;
	using Nd = std::shared_ptr<Node>;
	
	
	
	
	inline void AddAndIncrement(qi::symbols<char,int> * sym, std::string const& c, int * v)
	{
		sym->add(c,*v);
		(*v)++;
	}
	
	
	
	
	
	
	template<typename Iterator, typename Skipper = ascii::space_type>
	struct SystemParser : qi::grammar<Iterator, System(), Skipper>
	{
		
		
		SystemParser() : SystemParser::base_type(root_rule_)
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
			
			
			declarative_symbols.add("variable_group",0);
			declarative_symbols.add("hom_variable_group",1);
			declarative_symbols.add("variable",2);
			declarative_symbols.add("function",3);
			declarative_symbols.add("constant",4);
			declarative_symbols.add("parameter",5);
			declarative_symbols.add("implicit_parameter",6);
			declarative_symbols.add("pathvariable",7);
			declarative_symbols.add("random",8);
			declarative_symbols.add("random_real",9);
			
			
			root_rule_.name("root_rule_");
			root_rule_ =
			*(
				variable_group_ [phx::bind(&System::AddVariableGroup, _val, _1)]
				|
				hom_variable_group_ [phx::bind(&System::AddHomVariableGroup, _val, _1)]
				|
				variables_ [phx::bind(&System::AddUngroupedVariables, _val, _1)]
			)
			>>
			*(
				functions_ [phx::bind(&System::AddFunctions, _val, _1)]
			  )
			>>
			*(
				constants_ [phx::bind(&System::AddConstants, _val, _1)]
				  |
				parameters_ [phx::bind(&System::AddParameters, _val, _1)]
				  |
				implicit_parameters_ [phx::bind(&System::AddImplicitParameters, _val, _1)]
				  |
				path_variable_ [phx::bind(&System::AddPathVariable, _val, _1)]
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
			new_variable_ = unencountered_symbol_ [_val = make_shared_<Variable>() (_1)];
			
			
			
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
			unencountered_symbol_ = valid_variable_name_ - lexeme[( declarative_symbols | encountered_symbols )];
			// i am unsure about the use of lexeme in the above rule (unencountered_symbol).
			
			
			
			
			// get a string which fits the naming rules.
			valid_variable_name_.name("valid_variable_name_");
			valid_variable_name_ = +qi::alpha >> *(qi::alnum | qi::char_("[]_") );
			
			
			
			
			
			
			
			
			
//			debug(root_rule_);
//			debug(variable_group_);
//			debug(unencountered_symbol_);
//			debug(new_variable_);
//			debug(genericvargp_); (variable_group_) (valid_variable_name_)
//			BOOST_SPIRIT_DEBUG_NODES( (unencountered_symbol_) (new_variable_) (genericvargp_))
			
			using phx::val;
			using phx::construct;
			using namespace qi::labels;
			qi::on_error<qi::fail>
			( root_rule_ ,
			 std::cout<<
			 val("Error! Expecting ")<<
			 _4<<
			 val(" here: ")<<
			 construct<std::string>(_3,_2)<<
			 std::endl
			 );
			
			
		}
		
		
		
		
	private:
		
		// rule declarations.  these are member variables for the parser.
		qi::rule<Iterator, System(), Skipper > root_rule_;
		
		
		qi::rule<Iterator, std::vector<Var>(), Skipper > variable_group_, hom_variable_group_, variables_, implicit_parameters_;
		qi::rule<Iterator, std::vector<Var>(), Skipper > genericvargp_;
		
		qi::rule<Iterator, Var(), Skipper> path_variable_;
		qi::rule<Iterator, Var()> new_variable_;
		
		
		
		
		qi::rule<Iterator, std::vector<Fn>(), Skipper > functions_, constants_, parameters_;
		qi::rule<Iterator, std::vector<Fn>(), Skipper > genericfuncgp_;
		qi::rule<Iterator, Fn()>  new_function_;
		
		
		qi::rule<Iterator, std::string()> unencountered_symbol_;

		
		// the rule which determines valid variable names
		qi::rule<Iterator, std::string()> valid_variable_name_;
		
		
		// symbol declarations
		qi::symbols<char,Nd> encountered_symbols;
		qi::symbols<char,int> declarative_symbols;
		
		
		void MakeAndAddFunction(Fn & F, std::string str)
		{
			F = std::make_shared<Function>(str);
			encountered_symbols.add(str, F);
		}
		
	};
	
	
	
			
	
	
	template<typename Iterator>
	struct PartsParser : qi::grammar<Iterator, System(), boost::spirit::ascii::space_type>
	{
		
		
		
		PartsParser() : PartsParser::base_type(root_rule_) //,"SystemParser"
		{
			
			namespace phx = boost::phoenix;
			using qi::_1;
			using qi::_2;
			using qi::_3;
			using qi::_4;
			using qi::_val;
			using qi::eps;
			using qi::lit;
			
			
			
			
			
			
			
			
			
			///  declare variables to hold the function trees for defined symbols.
			std::vector<Nd> functions;
			std::vector<Nd> explicit_parameters;
			std::vector<Nd> subfunctions;
			std::vector<Nd> constants;
			std::vector<Var> variables;
			
			
			
			
			
			
			
			
			int variable_counter = 0;
			int path_variable_counter = 0;
			int function_counter = 0;
			int expl_para_counter = 0;
			int impl_para_counter = 0;
			int constant_counter = 0;
			int subfunction_counter = 0;
			
			
			
			qi::symbols<char,int> encountered_path_variable;
			
			
			
			qi::symbols<char,int> encountered_implicit_parameters;
			qi::symbols<char,int> encountered_explicit_parameters;
			
			
			qi::symbols<char,int> encountered_variables;
			encountered_variables.add("slartibartfast",0);  variable_counter++;
			
			qi::symbols<char,int> encountered_functions;
			
			qi::symbols<char,int> encountered_constant_functions;
			qi::symbols<char,int> encountered_subfunctions;
			
			
			
			
			VariableParser<Iterator> variable_parser;
			BrakeParser<Iterator> function_parser(&variables, encountered_variables);
			
			
			
			
			
			
			
			BOOST_SPIRIT_DEBUG_NODE(root_rule_);
			root_rule_.name("root_rule_");
			root_rule_ = eps >>
			*(vargp_ [phx::bind (
 						[](System & s, std::vector<Var> v )
 						{
 							std::cout << s.NumVariables() << std::endl;
 							s.AddVariableGroup(v);
 							std::cout << s.NumVariables() << std::endl;
 						}
 						,_val, _1) ]
			  );  // | definition_) >> ";";
			
			
			
			
			///////////////////
			//
			//  Declarations of things.
			//
			/////////////////////////
		
			
			
			
			
			
			BOOST_SPIRIT_DEBUG_NODE(vargp_);
			debug(vargp_);
			vargp_.name("vargp_");
			vargp_ %= "variable_group " > *( new_variable >> lit(',') )[phx::push_back(_val, _1)] > new_variable [phx::push_back(_val, _1)] > lit(';');
			
			
			
			BOOST_SPIRIT_DEBUG_NODE(new_variable);
			debug(new_variable);
			new_variable.name("new_variable");
			new_variable = unencountered_symbol_ [_val = make_shared_<bertini::Variable>()(_1)];
			
			
			// this rule gets a string.
			BOOST_SPIRIT_DEBUG_NODE(unencountered_symbol_);
			debug(unencountered_symbol_);
			unencountered_symbol_.name("unencountered_symbol");
			unencountered_symbol_ %= valid_variable_name_ - (encountered_variables);  // [_val = qi::_1]  [std::cout << _1 << std::endl]
			
			// |encountered_functions | encountered_explicit_parameters | encountered_implicit_parameters | encountered_constant_functions
			
			
			
			
			
			// get a string which fits the naming rules for Bertini, and Bertini2.
			BOOST_SPIRIT_DEBUG_NODE(valid_variable_name_);
			valid_variable_name_.name("valid_variable_name_");
			valid_variable_name_ %= +qi::alpha >> *(qi::alnum | qi::char_('_') | qi::char_('[') | qi::char_(']') );
			
			
			
			
			
			
			
			
//			variable_group_.name("variable_group_");
//			variable_group_ %=
//			"variable_group "
//			>>
//			unencountered_symbol_
//			 % ',' ;
			
			
//			[boost::phoenix::bind( [](qi::symbols<char,int> * sym, std::string c, int * v)
//								  {
//									  AddAndIncrement(sym, c, v);
//								  }, &encountered_variables, qi::_1, &variable_counter ) ]
			
			
//			implicit_parameters_.name("implicit_parameter_");
//			implicit_parameters_ = qi::lit("implicit_parameter ") >> variable_parser.valid_variable_name_
//			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
//								  {
//									  AddAndIncrement(sym, c, v);
//								  }, encountered_implicit_parameters, qi::_1, impl_para_counter ) ] % ',';
//			//[encountered_implicit_parameters.add(_1,impl_para_counter++)]
//			
//			
//			
//			
//			explicit_parameters_.name("explicit_parameter_");
//			explicit_parameters_ = qi::lit("explicit_parameter ") >> variable_parser.valid_variable_name_
//			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
//								  {
//									  AddAndIncrement(sym, c, v);
//								  }, encountered_explicit_parameters, qi::_1, expl_para_counter ) ] % ',' ;
//			
//			
//			constants_.name("constant_");
//			constants_ = (qi::lit("constant ") >> variable_parser.valid_variable_name_)
//			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
//								  {
//									  AddAndIncrement(sym, c, v);
//								  }, encountered_constant_functions, qi::_1, constant_counter ) ] % ',' ;
//			
//			
//			
//			functions_.name("function_");
//			functions_ = qi::lit("function ") >> variable_parser.valid_variable_name_
//			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
//								  {
//									  AddAndIncrement(sym, c, v);
//								  }, encountered_functions, qi::_1, function_counter ) ] % ',' ;
//			
//			
//			path_variable_.name("path_variable_");
//			path_variable_ =
//			qi::lit("path_variable ") >> variable_parser.valid_variable_name_ [_val = make_shared_<Variable>()(_1)] ;
			
			
			
			
			
		
			
			
			
			///////////////////
			//
			//  Definitions of things
			//
			////////////////////
			
			
			
//			definition_.name("definition_");
//			definition_ =
//			(subfunction_definition_ | function_definition_ )//| constant_definition_ | explicit_parameter_definition_)
//			>>
//			lit(";") ;
			
			
			
			
			
			
			
			
			
			
			
			
//			subfunction_definition_.name("subfunction_definition");
//			subfunction_definition_ = unencountered_symbol_ [ std::cout << qi::_1];
//			[ boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
//																				  {
//																					  AddAndIncrement(sym, c, v);
//																				  },
//																				  encountered_subfunctions, qi::_1, subfunction_counter ) ];
			
			
			
			
			
			// this actually makes a subfunction.
			
			
			
			
//			(unencountered_symbol_ >> lit("="));// >> function_parser);
//				[phx::bind( []
//													  (char found_symbol, Fn F, const qi::symbols<char, int> & sym, std::vector<Fn> & fns )
//													  {
//														  fns.push_back(F);
//														  auto bla = sym.find(found_symbol); // get a pointer to the object associated with the name of the function
//														  fns[*bla] = F;
//													  },
//													  _1,_2, phx::ref(encountered_subfunctions), phx::ref(subfunctions))];
			
			
			
			
//			
//			
//			
//			
//
			
//			
//			function_definition_.name("function_definition");
//			function_definition_ =
//			(encountered_functions >> lit("=") >> function_parser) [std::cout << "function " << _1 << " = " << _2 << std::endl];
			
			
//			[phx::bind(
//																			[]
//																			(const qi::symbols<char, int> & found_symbol, Nd F, const qi::symbols<char, int> & sym, std::vector<Nd> & fns )
//																			{
//																				fns.push_back(F); // this is so wrong
//																				auto bla = sym.find(found_symbol.name()); // get a pointer to the object associated with the name of the function
//																				fns[*bla] = F;
//																			},
//																			_1,_2, encountered_functions, functions)];
//
//			
//			
//			
//			constant_definition_.name("constant_definition");
//			constant_definition_ =
//			encountered_constant_functions >> lit("=") >> function_parser [phx::bind(
//																					 []
//																					 (const qi::symbols<char, int> & found_symbol, Fn F, const qi::symbols<char, int> & sym, std::vector<Fn> & fns )
//																					 {
//																						 fns.push_back(F);
//																						 auto bla = sym.find(found_symbol.name()); // get a pointer to the object associated with the name of the function
//																						 fns[*bla] = F;
//																					 },
//																					 _1,_2, encountered_constant_functions, constants)];
//			
//			
//			
//			
//			
//			
//			
//			explicit_parameter_definition_.name("explicit_parameter_definition");
//			explicit_parameter_definition_ =
//			encountered_explicit_parameters >> lit("=") >> function_parser [phx::bind(
//																					  []
//																					  (const qi::symbols<char, int> & found_symbol, Fn F, const qi::symbols<char, int> & sym, std::vector<Fn> & fns )
//																					  {
//																						  fns.push_back(F);
//																						  auto bla = sym.find(found_symbol.name()); // get a pointer to the object associated with the name of the function
//																						  fns[*bla] = F;
//																					  },
//																					  _1,_2, encountered_explicit_parameters, explicit_parameters)];
			
			
			
			
			
			
			

			
			
		}
		
		
		
		
		
		// rule declarations.  these are member variables for the parser.
		qi::rule<Iterator, System(), ascii::space_type > root_rule_;
		
		
		qi::rule<Iterator, ascii::space_type, std::vector<Var>()> vargp_;
		
		qi::rule<Iterator, ascii::space_type, Var()>  new_variable;
		
		qi::rule<Iterator, std::string(), ascii::space_type > unencountered_symbol_;// , ascii::space_type
		
		
		// the rule which determines valid variable names
		qi::rule<Iterator, std::string()> valid_variable_name_;
		
		
		qi::rule<Iterator, ascii::space_type, std::vector<Nd>() > declaration_;
		
		
		
		
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > explicit_parameters_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > implicit_parameters_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > path_variable_;
		qi::rule<Iterator, ascii::space_type > variable_group_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > constants_;
		qi::rule<Iterator, boost::spirit::unused_type > functions_;
		
		
		
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > definition_;
		
		
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > function_definition_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > constant_definition_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > explicit_parameter_definition_;
		qi::rule<Iterator, ascii::space_type, boost::spirit::unused_type > subfunction_definition_; //, boost::spirit::unused_type
		
		
		
		
		
		
		
		
		
	}; // struct SystemParser
	
	
	
	
	
}







#endif



