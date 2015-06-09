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

#include "system.hpp"
#include "function_tree/function_parsing.hpp"

namespace bertini {
	
	
	
	
	
	
	template<typename Iterator>
	struct SystemParser : qi::grammar<Iterator, System, boost::spirit::ascii::space_type>
	{
		
		// a few local using statements to reduce typing etc.
		using Fn = std::shared_ptr<Function>;
		using Var = std::shared_ptr<Variable>;
		
		
		SystemParser() : SystemParser::base_type(root_rule_) //,"SystemParser"
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
			std::vector<Fn> functions;
			std::vector<Fn> explicit_parameters;
			std::vector<Fn> subfunctions;
			std::vector<Fn> constants;
			std::vector<Var> variables;
			
			
			
			
			
			
			
			
			unsigned variable_counter = 0;
			unsigned path_variable_counter = 0;
			unsigned function_counter = 0;
			unsigned expl_para_counter = 0;
			unsigned impl_para_counter = 0;
			unsigned constant_counter = 0;
			unsigned subfunction_counter = 0;
			
			
			
			qi::symbols<char,int> encountered_path_variable;
			
			
			
			qi::symbols<char,int> encountered_implicit_parameters;
			qi::symbols<char,int> encountered_explicit_parameters;
			
			
			qi::symbols<char,int> encountered_variables;
			
			
			qi::symbols<char,int> encountered_functions;
			
			qi::symbols<char,int> encountered_constant_functions;
			qi::symbols<char,int> encountered_subfunctions;
			
			
			
			
			VariableParser<Iterator> variable_parser;
			BrakeParser<Iterator> function_parser(&variables, encountered_variables);
			
			
			
			
			
			
			
			
			root_rule_.name("root_rule_");
			root_rule_ = eps [_val = System()] >>
			*(declaration_ | definition_) ;
			
			
			
			///////////////////
			//
			//  Declarations of things.
			//
			/////////////////////////
			
			declaration_ =
			(implicit_parameters_ | path_variable_ | variable_group_ | functions_ | constants_ )// | hom_variable_group | variable)
			>>
			lit(";");
			
			
			
			
			implicit_parameters_.name("implicit_parameter_");
			qi::lit("implicit_parameter ") >> variable_parser.valid_variable_name_
			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
								  {
									  AddAndIncrement(sym, c, v);
								  }, encountered_implicit_parameters, qi::_1, impl_para_counter ) ] % ',';
			//[encountered_implicit_parameters.add(_1,impl_para_counter++)]
			
			
			
			
			explicit_parameters_.name("explicit_parameter_");
			qi::lit("explicit_parameter ") >> variable_parser.valid_variable_name_
			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
								  {
									  AddAndIncrement(sym, c, v);
								  }, encountered_explicit_parameters, qi::_1, expl_para_counter ) ] % ',' ;
			
			
			constants_.name("constant_");
			qi::lit("constant ") >> variable_parser.valid_variable_name_
			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
								  {
									  AddAndIncrement(sym, c, v);
								  }, encountered_constant_functions, qi::_1, constant_counter ) ] % ',' ;
			
			
			variable_group_.name("variable_group_");
			qi::lit("variable ") >> variable_parser.valid_variable_name_
			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
								  {
									  AddAndIncrement(sym, c, v);
								  }, encountered_variables, qi::_1, variable_counter ) ] % ',' ;
			
			
			functions_.name("function_");
			qi::lit("function ") >> variable_parser.valid_variable_name_
			[boost::phoenix::bind( [this](qi::symbols<char,int> & sym, std::string const& c, int & v)
								  {
									  AddAndIncrement(sym, c, v);
								  }, encountered_functions, qi::_1, function_counter ) ] % ',' ;
			
			
			path_variable_.name("path_variable_");
			path_variable_ =
			qi::lit("path_variable ") >> variable_parser.valid_variable_name_ [_val = make_shared_<Variable>()(_1)] ;
			
			
			
			
			
			
			
			
			
			
			
			
			
			///////////////////
			//
			//  Definitions of things
			//
			////////////////////
			
			
			
			
			definition_ %=
			(subfunction_definition_ | function_definition_ | constant_definition_ | explicit_parameter_definition_)
			>>
			lit(";") ;
			
			
			unencountered_symbol_.name("unencountered_symbol");
			unencountered_symbol_ = (
									 variable_parser.valid_variable_name_ -
									 (encountered_functions | encountered_explicit_parameters | encountered_implicit_parameters | encountered_constant_functions)
									 );
			
			
			
			
			subfunction_definition_.name("subfunction_definition");
//			subfunction_definition_ = // this actually makes a subfunction.
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
//			encountered_functions >> lit("=") >> function_parser [phx::bind(
//																			[]
//																			(const qi::symbols<char, int> & found_symbol, Fn F, const qi::symbols<char, int> & sym, std::vector<Fn> & fns )
//																			{
//																				fns.push_back(F);
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
			
			
			
			
			
			debug(root_rule_);
		}
		
		
		// rule declarations.  these are member variables for the parser.
		qi::rule<Iterator, System, ascii::space_type > root_rule_;
		
		
		
		qi::rule<Iterator, ascii::space_type > definition_;
		
		
		
		
		qi::rule<Iterator, ascii::space_type > explicit_parameters_;
		qi::rule<Iterator, ascii::space_type > implicit_parameters_;
		qi::rule<Iterator, ascii::space_type > path_variable_;
		qi::rule<Iterator, ascii::space_type > variable_group_;
		qi::rule<Iterator, ascii::space_type > constants_;
		qi::rule<Iterator, ascii::space_type > functions_;
		
		
		qi::rule<Iterator, ascii::space_type > declaration_;
		
		
		qi::rule<Iterator, ascii::space_type > function_definition_;
		qi::rule<Iterator, ascii::space_type > constant_definition_;
		qi::rule<Iterator, ascii::space_type > explicit_parameter_definition_;
		qi::rule<Iterator, ascii::space_type > subfunction_definition_;
		
		
		
		qi::rule<Iterator, std::string, ascii::space_type > unencountered_symbol_;
		
		
		
		
		void AddAndIncrement(qi::symbols<char,int> & sym, std::string const& c, int & v)
		{
			sym.add(c,v);
			v++;
		}
		
		
	}; // struct SystemParser
	
	
	
	
}







#endif



