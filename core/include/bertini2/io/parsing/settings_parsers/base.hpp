//This file is part of Bertini 2.
//
//bertini2/io/parsing/settings_parsers/base.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/settings_parsers/base.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/settings_parsers/base.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/settings_parsers/base.hpp
 
 \brief Provides the base type for config settings parsers
 */

#pragma once




#include "bertini2/io/parsing/qi_files.hpp"
#include "bertini2/io/parsing/number_rules.hpp"

namespace bertini {
	namespace parsing {
		namespace classic {

			namespace qi = ::boost::spirit::qi;
			namespace ascii = ::boost::spirit::ascii;
			
			
			/**
			 Qi Parser object for parsing config settings  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
			 This is a base case template, which we specialize elsewhere
			 */
			template<typename Iterator, typename Structure, typename Skipper = ascii::space_type> //boost::spirit::unused_type
			struct ConfigSettingParser : qi::grammar<Iterator, Structure(), Skipper>
			{ 
			};
			

		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini