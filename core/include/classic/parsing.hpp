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

			void SetInput(std::string new_input)
			{
				input_ = new_input;
			}

			void SetConfig(std::string new_config)
			{
				config_ = new_config;
			}

			void SetConfigInput(std::string c, std::string i)
			{
				config_ = c;
				input_ = i;
			}

			friend std::ostream& operator<<(std::ostream & out, SplitInputFile const& printme)
			{
				out << "--------config-----------\n\n" << printme.Config() << "\n\n-------input--------\n\n" << printme.Input();
				return out;
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

			\brief Qi Parser object for parsing text into the SplitInputFile class.  
			*/
			template<typename Iterator, typename Skipper = ascii::space_type> //boost::spirit::unused_type
			struct SplitFileInputConfig : qi::grammar<Iterator, SplitInputFile(), Skipper>
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
					using qi::char_;
					using qi::omit;
					using boost::spirit::lexeme;
					using boost::spirit::as_string;

					root_rule_.name("SplitFileInputConfig_root_rule");

					// "CONFIG\n\ntracktype:1;\n\nEND;\n\nINPUT\n\nvariable_group x, y;\nfunction f;\nf = x^2 + y^2 - 1;\n\nEND;\n\n";
					
					// (omit[lexeme[*(char_ - "CONFIG")]] > as_string[lexeme[*(char_ - "END;")]]  >> omit[*(char_ - "INPUT")] >> as_string[lexeme[*(char_ - "END;")]] >> omit[*char_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)];
					root_rule_ =  
								  (
								  	(
					              		(config_end_ > (input_end_ | end_ | input_ | no_decl_)) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
					              	)
					              	|
					              	((
					              	 	input_end_ | end_ | input_ | no_decl_
					              	) 
					              		[phx::bind(&SplitInputFile::SetInput, _val, _1)])
					              );
					//

					config_end_.name("config_end_");
					config_end_ = omit[*(char_ - "CONFIG")] >> ("CONFIG" > lexeme[*(char_ - "END;")] > "END;");

					input_end_.name("input_end_");
					input_end_ = omit[*(char_ - "INPUT")] >> "INPUT" >> end_;

					end_.name("end_");
					end_ = lexeme[*(char_ - "END;")] >> omit[*char_];

					input_.name("input_");
					input_ = omit[*(char_ - "INPUT")] >> "INPUT" >> *char_;

					no_decl_.name("no_decl_");
					no_decl_ = lexeme[*(char_)];


					// debug(root_rule_);debug(config_end_);debug(input_end_);debug(end_);debug(input_);debug(no_decl_);

					// BOOST_SPIRIT_DEBUG_NODES((root_rule_) 
					//                          (config_end_) 
					//                          (input_end_) (end_) (input_) (no_decl_) )



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
				}

			private:
				qi::rule<Iterator, SplitInputFile(), ascii::space_type > root_rule_;
				qi::rule<Iterator, ascii::space_type, std::string()> input_end_, end_, input_, no_decl_, config_end_;
			};

		} // re: namespace parsing

	} // re: namespace classic

} // re: namespace bertini






#endif

