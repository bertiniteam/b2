//This file is part of Bertini 2.
//
//bertini2/io/parsing/settings_parsers/random.hpp is free software: you can redistribute it
//and/or modify it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/io/parsing/settings_parsers/random.hpp is distributed in the hope that it will be
//useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/io/parsing/settings_parsers/random.hpp.  If not, see
//<http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
 \file bertini2/io/parsing/settings_parsers/random.hpp

 \brief Parser specialisation for RandomConfig — handles `randomseed: <uint>;`.
 */

#pragma once

#include "bertini2/io/parsing/settings_parsers/base.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"

namespace bertini {
	namespace parsing {
		namespace classic {

			template<typename Iterator, typename Skipper>
			/// \brief Parser for the RandomConfig settings block of classic Bertini input.
			struct ConfigSettingParser<Iterator, algorithm::RandomConfig, Skipper>
				: qi::grammar<Iterator, algorithm::RandomConfig(), Skipper>
			{
				ConfigSettingParser() : ConfigSettingParser::base_type(root_rule_, "config::Random")
				{
					namespace phx = boost::phoenix;
					using qi::_1;
					using qi::_val;
					using qi::lit;
					using qi::char_;
					using qi::ulong_;
					using boost::spirit::ascii::no_case;

					std::string seed_name = "randomseed";

					root_rule_.name("config::Random");
					root_rule_ = (random_seed_[phx::bind(
									[](algorithm::RandomConfig& S, unsigned long v) {
										S.random_seed = v;
									}, _val, _1)]
								  >> -no_setting_)
							   | no_setting_;

					all_names_ = no_case[seed_name] >> ':';

					random_seed_.name("random_seed_");
					random_seed_ = *(char_ - all_names_)
								 >> (no_case[seed_name] >> ':')
								 >> ulong_[_val = _1]
								 >> ';';

					no_setting_.name("no_setting_");
					no_setting_ = *(char_ - all_names_);
				}

			private:
				qi::rule<Iterator, algorithm::RandomConfig(), Skipper> root_rule_;
				qi::rule<Iterator, unsigned long(), Skipper> random_seed_;
				qi::rule<Iterator, Skipper, std::string()> no_setting_, all_names_;
			};

		} // namespace classic
	} // namespace parsing
} // namespace bertini
