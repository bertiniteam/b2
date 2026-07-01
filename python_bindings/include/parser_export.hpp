//This file is part of Bertini 2.
//
//python/include/parser_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/include/parser_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/include/parser_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
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
//  silviana amethyst
//  University of Wisconsin - Eau Claire
//  Spring 2018
//
//  python/include/parser_export.hpp:  Header file for exposing a method for parsing an input file to a system to python.



#pragma once
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
		// Route through the canonical C++ parse entry (bertini::parsing::classic::parse) rather
		// than duplicating the phrase_parse here: that one function owns the post-parse step that
		// emits the eager-bound functions into the System (EmitDeclaredFunctions).  Duplicating it
		// silently dropped every function (size-0 systems).
		template <typename ResultT, typename ParserT>
		ResultT Parser(std::string str)
		{
			// Treat input as UTF-8; drop a leading BOM so it is not parsed as a stray character.
			bertini::parsing::classic::StripUTF8BOM(str);
			ResultT res;
			bertini::parsing::classic::parse(str.begin(), str.end(), res);
			return res;
		};
		
		
		void ExportParsers();
	}
}

#endif
