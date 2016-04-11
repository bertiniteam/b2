//This file is part of Bertini 2.0.
//
// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//This file is part of Bertini 2.0.
//
// python/bertini_python.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/bertini_python.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/bertini_python.hpp.  If not, see <http://www.gnu.org/licenses/>.
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


#include <bertini2/function_tree/function_parsing.hpp>
#include <bertini2/system_parsing.hpp>



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
			def("parse_system", &Parser<System, SystemParser<std::string::const_iterator> >);
		};
		
	}
}

#endif
