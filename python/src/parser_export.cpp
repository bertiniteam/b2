//This file is part of Bertini 2.
//
//python/src/parser_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/src/parser_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/src/parser_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
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
//  Danielle Brake
//  UWEC
//  Spring 2018
//
//
//  python/src/parser_export.cpp



#include "parser_export.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini;
		

		
		void ExportParsers()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".parse");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("parse") = new_submodule;
			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Parsing functions, to turn strings into other things.";
			
			using namespace bertini::parsing;
			def("system", &Parser<System, classic::SystemParser<std::string::const_iterator> >);
		};
		
	}
}

