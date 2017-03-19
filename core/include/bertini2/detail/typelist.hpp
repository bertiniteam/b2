//This file is part of Bertini 2.
//
//typelist.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//typelist.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with typelist.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame


/**
\file detail/typelist.hpp

\brief Contains the detail/typelist helper type, for templated getting/setting/storage of config structs
*/


namespace bertini {
namespace detail {

/**
\brief A basic typelist, in the style of Modern C++ Design, but using variadic templates

To construct one, feed the types you want into the template arguments of the TypeList
*/
template <typename... Ts>
struct TypeList {
	using ToTuple = std::tuple<Ts...>;
	using ToTupleOfVec = std::tuple<Vec<Ts>...>;
	using ToTupleOfReal = std::tuple<typename Eigen::NumTraits<Ts>::Real...>;
};


}} // close namespaces





