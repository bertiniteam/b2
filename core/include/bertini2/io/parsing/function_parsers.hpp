//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/function_parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/function_parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/function_parsers.hpp
 
 \brief Provides the parsers for functions in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/function_rules.hpp"



namespace bertini {
	namespace parsing {
		namespace classic {
		
			
			
			template <typename Iterator>
			bool parse(Iterator first, Iterator last, std::shared_ptr<bertini::node::Node>& func)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				FunctionParser<Iterator> S;
				
				std::shared_ptr<bertini::node::Node>  f{};
				bool r = phrase_parse(first, last,
									  S,
									  space,
									  f);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				func = f;
				return r;
			}
		} // re: namespace classic
			
	}// re: namespace parsing
}// re: namespace bertini