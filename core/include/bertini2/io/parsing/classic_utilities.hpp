//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/classic_utilities.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/classic_utilities.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/classic_utilities.hpp
 
 \brief Provides utility functions for parsing classic bertini input files.
 */

#pragma once




#include "bertini2/io/parsing/qi_files.hpp"
#include "bertini2/io/file_utilities.hpp"



namespace bertini {
	namespace parsing {
		namespace classic {
			
			namespace qi = ::boost::spirit::qi;
			namespace ascii = ::boost::spirit::ascii;

			
			
			/**
			 \class SplitInputFile
			 
			 \brief Storage class for input file that splits up the config portion from the input portion.
			 
			 ## Use
			 
			 SplitFileInputConfig parser return an instance of this class.  The user can then parse the config and input portions separately using the appropriate parser.
			 
			 */
			class SplitInputFile
			{
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
				
				
				
				//            void StripComments()
				//            {
				//				std::string::const_iterator iter = config_.begin();
				//				std::string::const_iterator end = config_.end();
				//
				//				CommentStripper<std::string::const_iterator> S;
				//
				//				phrase_parse(iter, end, S,boost::spirit::ascii::space, config_);
				//
				//				iter = input_.begin();
				//				end = input_.end();
				//
				//				phrase_parse(iter, end, S,boost::spirit::ascii::space, input_);
				//            }
				
				
				friend std::ostream& operator<<(std::ostream & out, SplitInputFile const& printme)
				{
					out << "--------config-----------\n\n" << printme.Config() << "\n\n-------input--------\n\n" << printme.Input();
					return out;
				}
				
				
				
				
				
			private:
				std::string config_;
				std::string input_;
				
				bool readable_ = true; //Input file can be split accurately
				
			}; //re: SplitInputFile class

			
			
			
			/**
			 Qi Parser object for parsing text into the SplitInputfile class.  This ensures we can provide backwards compatibility with Bertini Classic input files.
			 
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
			 
			 CON END IN END  (1 1 1 1) 15
			 CON END IN --       (1 1 1 0) 14
			 CON END -- END   (1 1 0 1) 13
			 CON END -- --         (1 1 0 0)  12
			 CON -- IN END       (1 0 1 1) 11
			 CON -- IN  --           (1 0 1 0)  10
			 -- -- IN END            (0 0 1 1)  3
			 -- -- IN  --                 (0 0 1 0) 2
			 -- -- --END               (0 0 0 1) 1
			 --  --  --  --                 (0 0 0 0) 0
			 */
			template<typename Iterator, typename Skipper = ascii::space_type> //boost::spirit::unused_type
			struct SplitInputFileParser : qi::grammar<Iterator, SplitInputFile(), Skipper>
			{
				
				
				SplitInputFileParser() : SplitInputFileParser::base_type(root_rule_, "SplitInputFile")
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
					
					root_rule_.name("SplitInputFile_root_rule");
					
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
					config_ = no_case[*(char_ - "CONFIG") >> "CONFIG"];
					
					end_.name("end_");
					end_ = no_case[lexeme[*(char_ - "END;")] >> "END;"];
					
					input_.name("input_");
					input_ = no_case[lexeme[*(char_ - "INPUT")] >> "INPUT"];
					
					
					
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
					
					root_rule_ = eps[_val = ""] >> *line_[_val = _val + _1 + "\n"] >> -last_line_[_val = _val + _1 + "\n"];//+line_ | qi::eoi;
					
					
					line_.name("line_of_commented_input");
					line_ = lexeme[*(char_ - eol - "%") >> omit[-( "%" >> lexeme[*(char_ - eol)] )] >> (eol ) ];
					
					last_line_.name("line_of_commented_input_with_no_eol");
					last_line_ = lexeme[*(char_ - "%") >> omit[-( "%" >> lexeme[*(char_ - eol)] )]] ;
					
					
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
			
			
			
			SplitInputFile ParseInputFile(std::string str_input_file)
			{
				std::string::const_iterator iter = str_input_file.begin();
				std::string::const_iterator end = str_input_file.end();
				
				std::string uncommented_str_input;
				CommentStripper<std::string::const_iterator> commentparser;
				phrase_parse(iter, end, commentparser,boost::spirit::ascii::space, uncommented_str_input);
				
				iter = uncommented_str_input.begin();
				end = uncommented_str_input.end();
				
				SplitInputFileParser<std::string::const_iterator> spliter;
				
				SplitInputFile input_file;
				phrase_parse(iter, end, spliter,boost::spirit::ascii::space, input_file);
				
				return input_file;
			}
			


			/**
			 \brief Function for splitting a Bertini Classic style input file into `config` and `input`.
			 */
			std::tuple<std::string, std::string> SplitIntoConfigAndInput(Path const& input_file)
			{
				auto file_as_string = FileToString(input_file);
				
				SplitInputFileParser<std::string::const_iterator> parser;
				SplitInputFile config_and_input;
				std::string::const_iterator iter = file_as_string.begin();
				std::string::const_iterator end = file_as_string.end();
				phrase_parse(iter, end, parser, boost::spirit::ascii::space, config_and_input);

				return std::make_tuple(config_and_input.Config(), config_and_input.Input());
			}


			/**
			 \brief Function for splitting a Bertini Classic style input file into `config` and `input`.
			 */
			void SplitIntoConfigAndInput(std::string & config_section, std::string & input_section, Path const& input_file)
			{
				std::tie(config_section, input_section) = SplitIntoConfigAndInput(input_file);
			}


		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini