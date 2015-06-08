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
		
		
		SystemParser() : SystemParser::base_type(root_rule_,"SystemParser")
		{
			namespace phx = boost::phoenix;
			using qi::_1;
			using qi::_2;
			using qi::_3;
			using qi::_4;
			using qi::_val;
			using qi::eps;
			using qi::lit;
			
			
			root_rule_.name("root_rule_");
			root_rule_ = eps [_val = System()];
			
			
			
			independent_symbol_.name("independent_symbol_");
			
			
			
			implicit_parameter_.name("implicit_parameter_");
			
			
			path_variable_.name("path_variable_");
			
			
			
			
			
			
			
			dependent_symbol_.name("dependent_symbol_");
			
			
			
			function_.name("function");
			
			
			subfunction_.name("subfunction");
			
			
			explicit_parameter_.name("explicit_parameter");
			
			
			
			
			
			debug(root_rule_);
			debug(independent_symbol_);
			debug(dependent_symbol_);
		}
		
	
		// rule declarations.  these are member variables for the parser.
		qi::rule<Iterator, System, ascii::space_type > root_rule_;
		
		
	
		
		qi::rule<Iterator, std::shared_ptr<Variable>, ascii::space_type > independent_symbol_;
		qi::rule<Iterator, std::shared_ptr<Variable>, ascii::space_type > implicit_parameter_;
		qi::rule<Iterator, std::shared_ptr<Variable>, ascii::space_type > path_variable_;
		
		qi::rule<Iterator, std::shared_ptr<Function>, ascii::space_type > dependent_symbol_;
		qi::rule<Iterator, std::shared_ptr<Function>, ascii::space_type > function_;
		qi::rule<Iterator, std::shared_ptr<Function>, ascii::space_type > subfunction_;
		qi::rule<Iterator, std::shared_ptr<Function>, ascii::space_type > explicit_parameter_;
		
		
		
	}; // struct SystemParser
	
	
	
	
}







#endif



