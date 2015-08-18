//This file is part of Bertini 2.0.
//
//parsing.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//parsing.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with parsing.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
// parsing.hpp: provides some parsers for classic mode input files.


#ifndef BERTINI_CLASSIC_PARSING
#define BERTINI_CLASSIC_PARSING

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/include/adapted.hpp>

#define BOOST_SPIRIT_USE_PHOENIX_V3 1


#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>


#include <iostream>
#include <string>


namespace bertini
{
	namespace classic
	{
		class SplitInputFile
		{
			std::string config_;
			std::string input_;
		public:
			
			std::string Config() const
			{
				return config_;
			}

			std::string Input() const
			{
				return input_;
			}

		};
		namespace parsing
		{
			/**
			Qi Parser object for parsing text into the System class.  This ensures we can provide backwards compatibility with Bertini Classic input files.

			To use this parser, construct an object of its type, then use it to parse.

			\code
			System sys;
			std::string str = "variable_group x, y, z; \nfunction f1, f2;\n  f1 = x*y*z;\n f2 = x+y+z;\n";

			std::string::const_iterator iter = str.begin();
			std::string::const_iterator end = str.end();


			bertini::SystemParser<std::string::const_iterator> S;


			bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);

			\endcode

			\brief Qi Parser object for parsing text into the System class.  
			*/
			template<typename Iterator, typename Skipper = ascii::space_type>
			struct SplitFileInputConfig : qi::grammar<Iterator, SplitInputFile, Skipper>
			{
				
				
				SplitFileInputConfig() : SplitFileInputConfig::base_type(root_rule_, "SplitFileInputConfig")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_2;
					using qi::_3;
					using qi::_4;
					using qi::_val;
					using qi::eps;
					using qi::lit;
					using boost::spirit::lexeme;


				}

			private:
				qi::rule<Iterator, SplitInputFile, ascii::space_type > root_rule_;
			};

		} // re: namespace parsing

	} // re: namespace classic

} // re: namespace bertini






#endif

