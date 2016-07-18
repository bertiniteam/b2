//This file is part of Bertini 2.
//
//amp_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// James Collins, West Texas A&M University


/**
 \file config_parsing.hpp
 
 \brief Contains qi parsers and functions needed to parse the config settings from a classic Bertini input file
 */

#ifndef config_parsing_hpp
#define config_parsing_hpp

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>



#define BOOST_SPIRIT_USE_PHOENIX_V3 1


#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi_no_case.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include <bertini2/Tracking/tracking_config.hpp>


#include <iostream>
#include <string>


namespace bertini
{
	namespace classic
	{
		namespace parsing
		{
			using namespace bertini::tracking;
			
			
			/**
			 Qi Parser object for parsing config settings  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
			 To use this parser, construct an object of its type, then use it to parse.
			 
			 \code
			 System sys;
			 std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";
			 
			 std::string::const_iterator iter = str.begin();
			 std::string::const_iterator end = str.end();
			 
			 
			 bertini::SystemParser<std::string::const_iterator> S;
			 
			 
			 bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);
			 
			 \endcode
			 
			 \brief Qi Parser object for parsing config file to determine MPType setting.
			 
			 */
			template<typename Iterator, typename Skipper = ascii::space_type> //boost::spirit::unused_type
			struct ConfigPrecisionTypeParser : qi::grammar<Iterator, PrecisionType(), Skipper>
			{
				
				
				ConfigPrecisionTypeParser() : ConfigPrecisionTypeParser::base_type(root_rule_, "ConfigPrecisionType")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;
					using boost::spirit::ascii::no_case;
					
					precisiontype_.add("0", PrecisionType::Fixed);
					precisiontype_.add("1", PrecisionType::Adaptive);
					precisiontype_.add("2", PrecisionType::Adaptive);


	
					
					root_rule_.name("ConfigPrecisionType_root_rule");
					
					root_rule_ = -(omit[config_name_] > precisiontype_ > ';');

					config_name_.name("precisiontype_");
					config_name_ = no_case[lexeme[*(char_ - "mptype:")] >> "mptype:"];

					

					
					using phx::val;
					using phx::construct;
					using namespace qi::labels;
					qi::on_error<qi::fail>
					( root_rule_ ,
					 std::cout<<
					 val("config/input split parser could not complete parsing. Expecting ")<<
					 _4<<
					 val(" here: ")<<
					 construct<std::string>(_3,_2)<<
					 std::endl
					 );
					
//					debug(root_rule_);
					//
					//
					//
					//					debug(config_name_);

	
				}
				
				
			private:
				qi::rule<Iterator, PrecisionType(), ascii::space_type > root_rule_;
				qi::rule<Iterator, ascii::space_type, std::string()> config_name_;
				
				qi::symbols<char,PrecisionType> precisiontype_;


//				struct precisiontype : qi::symbols<unsigned, std::string>
//				{
//					precisiontype()
//					{
//						add(1,"h");
//					}
//				}precisiontype_;
			};
			
			
		} // re: namespace parsing
		
	} // re: namespace classic
	
} // re: namespace bertini


#endif
