//This file is part of Bertini 2.
//
//bertini2/blackbox/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017-2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/blackbox/config.hpp 

\brief Configs for the blackbox whatnot
*/

#pragma once

namespace bertini{
	namespace blackbox{


namespace type{
enum class Start{ TotalDegree, MHom, User};
enum class Tracker{ FixedDouble, FixedMultiple, Adaptive};
enum class Endgame{ PowerSeries, Cauchy};
}


template<typename T>
	struct StorageSelector{};

template<>
	struct StorageSelector<start_system::User>
	{
		using ShouldClone = typename std::false_type;
		// typename algorithm::StorageSelector<StartType>::Storage
		// using Storage = typename policy::RefToGiven<System, start_system::User>;
	};

template<>
	struct StorageSelector<start_system::TotalDegree>
	{
		// using Storage = typename policy::CloneGiven<System, start_system::TotalDegree>;
		// typedef policy::CloneGiven Storage;
		using ShouldClone = typename std::true_type;
		// using Storage = typename policy::CloneGiven;
	};

template<>
	struct StorageSelector<start_system::MHomogeneous>
	{
		using ShouldClone = typename std::true_type;
		// template < typename T, typename S>
		// using Storage = typename policy::CloneGiven<T,S>;
	};


	}
}
