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



#include "bertini2/detail/typelist.hpp"

#include "bertini2/io/parsing/settings_parsers/tracking.hpp"
#include "bertini2/io/parsing/settings_parsers/endgames.hpp"
#include "bertini2/io/parsing/settings_parsers/algorithm.hpp"




namespace bertini {
	namespace parsing {
		
		namespace classic {
			
			/**
			 \brief Helper function to fill a single configuration struct by parsing a config input file.
			 
			 \param config_str The comment-stripped configuration string from a Bertini classic input file.
			 
			 \tparam ConfigT The config ConfigT type
			 \tparam RT Real number type
			 
			 \returns The config struct filled with data from the input file.
			 */
			
			template<typename ConfigT>
			ConfigT FillConfigStruct(std::string const& config_str)
			{
				std::string::const_iterator iter = config_str.begin();
				std::string::const_iterator end = config_str.end();
				ConfigT settings;
				ConfigSettingParser<std::string::const_iterator, ConfigT> parser;
				auto parse_success = phrase_parse(iter, end, parser,boost::spirit::ascii::space, settings);
				if (!parse_success || iter!=end)
					throw std::runtime_error("failed to parse into config struct from file");

				return settings;
			}




			/**
			\brief Base variadic parser for parsing many config structs at once.

			A specialization for a typelist of configs appears below.

			\tparam RT Real number type
			\tparam Ts Configuration structures to be filled by the parser
			*/
			template<typename ...Ts>
			struct ConfigParser
			{
				/**
				 The primary idea for filling a tuple using variadic templates comes from:
				 http://stackoverflow.com/questions/10014713/build-tuple-using-variadic-templates
				 
				 
				 \brief Reads in a comment-stripped, config portion of a Bertini classic input file.  Parses the config settings and returns the structures passed into the template parameters with the relevant config settings.
				 
				 \param config The string containing the comment-stripped config portion of the Bertini classic input file.
				 
				 \return A tuple containing all the required config structures.
				 */
				static
				std::tuple<Ts...> Parse(std::string const& config)
				{
					return std::make_tuple<Ts...>(FillConfigStruct<Ts>(config)...);
				}

			};


			/**
			\brief Specialization of ConfigParser for a single config struct

			\tparam ConfigT The config ConfigT type
			\tparam RT Real number type
			*/
			template<typename ConfigT>
			struct ConfigParser <ConfigT>
			{	

				/**
				 \brief Fill a single configuration struct by parsing a config input file.
				 
				 \param config The comment-stripped configuration string from a Bertini classic input file.
				 
				 \returns The config struct filled with data from the input file.
				 */
				static
				ConfigT Parse(std::string const& config)
				{
					return FillConfigStruct<ConfigT>(config);
				}
			};


			/**
			\brief Specialization of ConfigParser for a typelist, returning a thing passed down from the base variadic case.

			\tparam ConfigT The config ConfigT type
			\tparam RT Real number type
			*/
			template<typename ...Ts>
			struct ConfigParser<detail::TypeList<Ts...>>
			{	
				static
				auto Parse(std::string const& config)
				{
					return ConfigParser<Ts...>::Parse(config);
				}
			};
		} // re: namespace classic
		
	}// re: namespace parsing
}// re: namespace bertini