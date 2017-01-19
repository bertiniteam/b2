//This file is part of Bertini 2.
//
//configured.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//configured.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with configured.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Daniel Brake
// University of Notre Dame
//

//  detail/configured.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file detail/configured.hpp

\brief Contains the detail/configured helper type, for templated getting/setting/storage of config structs
*/

#pragma once

#include "bertini2/detail/is_template_parameter.hpp"

namespace bertini {

	namespace detail{


		/**
		This templated struct allows one to provide a variable typelist of config structs (or anything else for that matter), and provides a way to get and set these configs by looking up the structs by type.

		I will add lookup by index if it is needed.  If you need to look up by index, try type first.  If your types are non-unique, consider making a struct to hold the things you are storing.
		*/
		template<typename ...Ts>
		struct Configured
		{
			std::tuple<Ts...> configuration_;

			template<typename ... T>
			Configured(T const& ...t) : configuration_(t...)
			{}


			template<typename T, typename = typename std::enable_if<IsTemplateParameter<T,Ts...>::value>::type>
			const T& Get() const
			{
				return std::get<T>(configuration_);
			}

			template<typename T, typename = typename std::enable_if<IsTemplateParameter<T,Ts...>::value>::type>
			void Set(T const& t)
			{
				std::get<T>(configuration_) = t;
			}


		}; // Configured

		#define FORWARD_GET_CONFIGURED \
		template <typename T> \
		const T& Get() const \
		{ return Config::template Get<T>();}

	} // namespace detail
}