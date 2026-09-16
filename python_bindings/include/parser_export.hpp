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


#include <type_traits>

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
			// accept a full Bertini 1 classic file, which is exactly what
			// System.to_classic_input() emits -- see #396
			bertini::parsing::classic::StripClassicFileWrappers(str);
			ResultT res;
			// The return value is NOT optional to check.  parse() answers false when the
			// grammar matched nothing or did not consume the whole input, and it leaves the
			// result untouched -- so ignoring it handed Python a structurally valid, entirely
			// EMPTY System and reported nothing.  A caller round-tripping through text then
			// carried on with a system of zero functions.  See #396.
			if (!bertini::parsing::classic::parse(str.begin(), str.end(), res))
			{
				std::string shown = str.size() > 80 ? str.substr(0, 80) + "..." : str;
				throw std::runtime_error(
					"[SystemParser] could not parse a system from the input.  Expected classic "
					"declarations (variable_group / function / definitions), optionally wrapped "
					"in a Bertini 1 `INPUT ... END;` section with an optional leading "
					"`CONFIG ... END;` section.  Input began: \"" + shown + "\"");
			}
			// parse() also answers TRUE for input that declares NOTHING -- an empty string, or
			// an `INPUT ... END;` section with an empty body -- because matching zero
			// declarations is a successful match of the grammar.  The caller still ends up
			// holding an empty System they did not ask for, which is the whole complaint in
			// #396, so refuse that too.  A system with variables but no functions is left
			// alone: that is a real, if unusual, declaration.
			if constexpr (std::is_same<ResultT, System>::value)
			{
				if (res.NumVariables()==0 && res.NumTotalFunctions()==0)
					throw std::runtime_error(
						"[SystemParser] the input declared no variables and no functions, so "
						"there is no system to return.  Expected classic declarations "
						"(variable_group / function / definitions), optionally wrapped in a "
						"Bertini 1 `INPUT ... END;` section.");
			}
			return res;
		};
		
		
		void ExportParsers();
	}
}

#endif
