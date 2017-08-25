//This file is part of Bertini 2.
//
//is_template_parameter.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//is_template_parameter.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with is_template_parameter.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame
//

//  detail/is_template_parameter.hpp
//
//  copyright 2016
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring 2016

/**
\file detail/is_template_parameter.hpp

\brief Contains the detail/is_template_parameter helper, for checking whether a type is a template parameter in a variable template
*/

#pragma once
#include "bertini2/detail/typelist.hpp"

namespace bertini {

	namespace detail{

	/**
	Base template, which says 'no', it's not in the pack.
	*/
	template <typename...>
	struct IsTemplateParameter {
	    static constexpr bool value = false;
	};

	/**
	Specialized template, which says 'yes', if it is in the pack, by recursively expanding the pack.
	*/
	template <typename F, typename S, typename... T>
	struct IsTemplateParameter<F, S, T...> {
	    static constexpr bool value =
	        std::is_same<F, S>::value || IsTemplateParameter<F, T...>::value;
	        // either it is the same as the first one in the pack, or we need to expand to the right in the pack.
	};


	template < typename T, typename ...Ts>
	struct IsTemplateParameter<T, TypeList<Ts...>>
	{
		static constexpr bool value = IsTemplateParameter<T, Ts...>::value;
	};



	} // namespace detail


}