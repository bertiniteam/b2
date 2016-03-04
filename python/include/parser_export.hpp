//
//  parser_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/10/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_parser_export_hpp
#define Xcode_b2_parser_export_hpp
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
