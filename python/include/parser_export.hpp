//This file is part of Bertini 2.
//
//python/parser_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/parser_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/parser_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/parser_export.hpp:  Header file for exposing a method for parsing an input file to a system to python.




#ifndef BERTINI_PYTHON_PARSER_EXPORT_HPP
#define BERTINI_PYTHON_PARSER_EXPORT_HPP
#include <boost/spirit/include/qi.hpp>


#include <bertini2/io/parsing/function_parsers.hpp>
#include <bertini2/io/parsing/system_parsers.hpp>


#include "python_common.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini;
		
		/////////////  Parser Exposure  /////////////////////
		template <typename ResultT, typename ParserT>
		ResultT Parser(std::string str)
		{
			ResultT res;
			std::string::const_iterator iter = str.begin();
			std::string::const_iterator end = str.end();
			ParserT P;
			phrase_parse(iter, end, P, boost::spirit::ascii::space, res);
			
			return res;
		};
		
		
		
		
		void ExportParsers()
		{
			using namespace bertini::parsing;
			def("parse_system", &Parser<System, classic::SystemParser<std::string::const_iterator> >);
		};
		
	}
}

#endif
