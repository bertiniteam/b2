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
			
			
			
			
			
			
			
		}
		
	
		// rule declarations.  these are member variables for the parser.
		qi::rule<Iterator, System, ascii::space_type > root_rule_;
		
		
	}; // struct SystemParser
	
	
	
	
}







#endif



