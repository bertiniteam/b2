//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/system_parsers.hpp
 
 \brief Provides the parsers for systems in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/system_rules.hpp"



namespace bertini {
	namespace parsing {
		
		namespace classic {
		
		
		
			template <typename Iterator>
			static bool parse(Iterator first, Iterator last, System& sys)
			{
				using boost::spirit::qi::double_;
				using boost::spirit::qi::_1;
				using boost::spirit::qi::phrase_parse;
				using boost::spirit::ascii::space;
				using boost::phoenix::ref;
				
				SystemParser<Iterator> S;
				
				std::shared_ptr<Node>  s{};
				bool r = phrase_parse(first, last,
									  S,
									  space,
									  s);
				
				if (!r || first != last) // fail if we did not get a full match
					return false;
				
				std::swap(s,sys);
				return r;
			}
		
		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini