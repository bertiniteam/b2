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
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// configIni_parser.hpp: Converts configuration settings into ini format


#ifndef Function_Tree_Test_configIni_parse_hpp
#define Function_Tree_Test_configIni_parse_hpp


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
    namespace settings
    {
        namespace parsing
        {
            /**
             Qi Parser object for converting Bertini configuration settings into ini format so they can be read by Boost.ProgramOptions.
             
             To use this parser, construct an object of its type, then use it to parse.
             
             \code
             System sys;
             std::string str = "TrackType: 1\n MPType: 2";
             
             std::string::const_iterator iter = str.begin();
             std::string::const_iterator end = str.end();
             
             
             bertini::SystemParser<std::string::const_iterator> S;
             
             
             bool s = phrase_parse(iter, end, S,boost::spirit::ascii::space, sys);
             
             \endcode
             
             \brief Qi Parser object converting settings to ini format.
             */
            
            template<typename Iterator, typename Skipper = ascii::space_type> //boost::spirit::unused_type
            struct ConfigToIni : qi::grammar<Iterator, std::string(), Skipper>
            {
                
                
                ConfigToIni() : ConfigToIni::base_type(root_rule_, "ConfigToIni")
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
                    
                    root_rule_.name("ConfigToIni_root_rule");
                    
                    root_rule_ = eps[_val = ""] >> *line_[_val = _val + _1+"\n"] >> -last_line_[_val = _val + _1];
                    
                    
                    line_.name("line_of_settings_input");
                    line_ = (eps >> colon_ >> yes_eol_)[_val = _1 + "=" + _2];
                    
                    last_line_.name("line_of_settings_input_with_no_eol");
                    last_line_ = (eps >> colon_ >> no_eol_)[_val = _1 + "=" + _2];
                    
                    colon_.name("Set of characters with no : or eol");
                    colon_ = *(char_ -  ":") >> omit[":"];
                    
                    yes_eol_.name("Set of characters with eol");
                    yes_eol_ = lexeme[*(char_ - eol) >> eol];

                    no_eol_.name("Set of characters with no eol");
                    no_eol_ = lexeme[*(char_)];
                    
                    //                     debug(root_rule_);
                    //                     debug(line_);
                    //
                    //                     BOOST_SPIRIT_DEBUG_NODES((root_rule_)
                    //                                              (line_) )
                    
                    
                    
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
                qi::rule<Iterator, ascii::space_type, std::string()> line_, last_line_, colon_, yes_eol_, no_eol_;
                
            }; //end struct ConfigToIni
            
        } // end parsing namespace
    }// end settings namespace
}// end bertini namespace






#endif
