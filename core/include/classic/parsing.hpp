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
            
            bool readable_ = true; //Input file can be split accurately
		public:

			std::string Config() const
			{
				return config_;
			}

			std::string Input() const
			{
				return input_;
			}
            
            bool Readable() const
            {
                return readable_;
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
            
            void SetReadable(bool read)
            {
                readable_ = read;
            }
            
            
            
            void StripComments()
            {
                
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

					root_rule_ = both_ | unreadable_ |only_input_ ;

					// NOTE ON CONVENTIONS CURRENTLY USED:
					// the various rules are constructed to obtain the text BEFORE the named marker, so e.g, rule config_ gets all text in the string before 'CONFIG'.  similarly for 'END;' and 'INPUT'.

					both_.name("have_both_config_and_input");
					both_ = (
		              		 (omit[config_] >> end_ >> omit[input_] >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]// 15
		              		 |
		              		 (omit[config_] >> end_ >> omit[input_] >> no_decl_) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)] // 14
		              		 |
		              		 (omit[config_] >> end_ >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)] // 13
		              		 |
		              		 (omit[config_] >> input_ >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)] // 11
		              		 |
		              		 (omit[config_] >> end_ >> no_decl_) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)] // 12
		              		 |
		              		 (omit[config_] >> input_ >> no_decl_) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)] // 10
		              		)
		              		;

		              		// the rule for number 11 does not currently correctly parse text which matches it...  why???

		              		//this attempt for the both_ rule doesn't work because the ends are optional...
		              		// (omit[config_] >> end_ >> omit[input_] >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
		              		//  |
		              		//  (omit[config_] >> omit[!end_] >> input_ >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
		              		//  |
		              		//  (omit[config_] >> end_ >> omit[!input_] >> end_ >> omit[no_decl_]) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
		              		//  |
		              		//  (omit[config_] >> omit[!end_] >> input_ >> no_decl_) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
		              		//  |
		              		//  (omit[config_] >> end_ >> omit[!input_] >> omit[!end_] >> no_decl_) [phx::bind(&SplitInputFile::SetConfigInput, _val, _1, _2)]
		              		// )


		            only_input_.name("have_only_input");
					only_input_ = (
				              	 		(omit[input_] >> end_ >> omit[no_decl_]) // 3
				              	 		|
				              	 		(omit[input_] >> no_decl_) // 2
				              	 		|
				              	 		(end_ >> omit[no_decl_])  // 1
				              	 		|
				              	 		(no_decl_) // 0
				              		) 
				              		[phx::bind(&SplitInputFile::SetInput, _val, _1)];

                    unreadable_.name("unreadable_input"); //9 and 4 same as 12 and 1 to parser
                    unreadable_ = (
                                   (omit[config_] >> no_decl_) // 8
                                   |
                                   (end_ >> omit[input_] >> end_ >> omit[no_decl_])  // 7
                                   |
                                   (end_ >> omit[input_]  >> omit[no_decl_])  // 6
                                   |
                                   (end_ >>  end_ >> omit[no_decl_])  // 5
                                   )
                                    [phx::bind(&SplitInputFile::SetReadable, _val, false)];

					config_.name("config_");
					config_ = *(char_ - "CONFIG") >> "CONFIG";

					end_.name("end_");
					end_ = lexeme[*(char_ - "END;")] >> "END;";

					input_.name("input_");
					input_ = lexeme[*(char_ - "INPUT")] >> "INPUT";



					no_decl_.name("no_decl_");
					no_decl_ = lexeme[*(char_)];


					// debug(root_rule_);
					// debug(both_); debug(only_input_);
					// debug(config_);
					// debug(end_);
					// debug(input_);
					// debug(no_decl_);

					// BOOST_SPIRIT_DEBUG_NODES((root_rule_) 
					//                          (both_) (only_input_)
					//                          (config_) (end_) (input_) 
					//                          (no_decl_) )



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
				qi::rule<Iterator, SplitInputFile(), ascii::space_type > root_rule_, both_, only_input_, unreadable_;
				qi::rule<Iterator, ascii::space_type, std::string()> end_, input_, no_decl_;
				qi::rule<Iterator, ascii::space_type> config_;
			};

            
            
            
            
            
            /**
             Qi Parser object for removing comments from a line of text and then returns the uncommented line.
             
             To use this parser, construct an object of its type, then use it to parse.
             
             */
            template<typename Iterator, typename Skipper = ascii::space_type> //boost::spirit::unused_type
            struct CommentStripper : qi::grammar<Iterator, std::string(), Skipper>
            {
                
                
                CommentStripper() : CommentStripper::base_type(root_rule_, "CommentStripper")
                {
                    namespace phx = boost::phoenix;
                    using qi::_1;
                    using qi::_2;
                    using qi::_3;
                    using qi::_4;
                    using qi::_val;
                    using qi::eol;
                    using qi::eoi;
                    using qi::eps;
                    using qi::lit;
                    using qi::char_;
                    using qi::omit;
                    using boost::spirit::lexeme;
                    using boost::spirit::as_string;
                    
                    root_rule_.name("CommentStripper_root_rule");
                    
                    root_rule_ = eps[_val = ""] >> *line_[_val = _val + _1 + "=" + _2 + "\n"] >> -last_line_[_val = _val + _1 + "\n"];//+line_ | qi::eoi;
                    
                   
                    line_.name("line_of_commented_input");
                    line_ = lexeme[*(char_ - eol - "%") >> omit[-( "%" >> lexeme[*(char_ - eol)] )] >> (eol ) ];
                    
                    last_line_.name("line_of_commented_input_with_no_eol");
                    last_line_ = lexeme[*(char_ - eol - ":") >> omit[":"] >> lexeme[*(char_ - eol)] ] ;
                    
                    
//                     debug(root_rule_);
//                     debug(line_);
//                     debug(last_line_);
//
//                     BOOST_SPIRIT_DEBUG_NODES((root_rule_) 
//                                              (line_) (last_line_) )
                    
                    
                    
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
                qi::rule<Iterator, std::string(), ascii::space_type > root_rule_;
                qi::rule<Iterator, ascii::space_type, std::string()> line_, last_line_;
            };

		} // re: namespace parsing

	} // re: namespace classic

} // re: namespace bertini






#endif

