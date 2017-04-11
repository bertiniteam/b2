//This file is part of Bertini 2.
//
//bertini2/io/parsers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/settings_parsers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/settings_parsers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
 \file bertini2/io/parsing/settings_parsers.hpp
 
 \brief Provides the parsers for settings in bertini2.
 */

#pragma once




#include "bertini2/io/parsing/settings_rules.hpp"



namespace bertini {
	namespace parsing {
		
		namespace classic {
			
			
			/**
			 \brief Helper function to fill a single configuration struct by parsing a config input file.
			 
			 \param config_str The comment-stripped configuration string from a Bertini classic input file.
			 
			 \tparam Structure The config structure type
			 \tparam RT Real number type
			 
			 \returns The config struct filled with data from the input file.
			 */
			
			template<typename Structure, typename RT>
			Structure FillConfigStruct(std::string config_str)
			{
				std::string::const_iterator iter = config_str.begin();
				std::string::const_iterator end = config_str.end();
				Structure settings;
				ConfigSettingParser<std::string::const_iterator, Structure, RT> parser;
				// TODO: error handling if parsing goes wrong
				bool parsed = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
				
				return settings;
			}
			
			
			
			
			/**
			 This idea for filling a tuple using variadic templates comes from:
			 http://stackoverflow.com/questions/10014713/build-tuple-using-variadic-templates
			 
			 
			 \brief Reads in a comment-stripped, config portion of a Bertini classic input file.  Parses the config settings and returns the structures passed into the template parameters with the relevant config settings.
			 
			 \param config_str The string containing the comment-stripped config portion of the Bertini classic input file.
			 
			 \tparam RT Real number type
			 \tparam Structs A configuration structure to be filled by the parser
			 
			 \return A tuple containing all the required config structures.
			 */
			template<typename RT, typename... Structs>
			std::tuple<Structs...> GetConfigSettings(std::string config_str)
			{
				return std::make_tuple<Structs...>(FillConfigStruct<Structs,RT>(config_str)...);
			}
			
		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini