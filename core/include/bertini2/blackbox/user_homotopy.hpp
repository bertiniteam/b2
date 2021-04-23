//This file is part of Bertini 2.
//
//bertini2/blackbox/switches_zerodim.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/switches_zerodim.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/switches_zerodim.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017-2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/blackbox/switches_zerodim.hpp 

\brief A sequence of switches for getting a particular instantiation of a ZeroDim algorithm, based on runtime options.
*/



#pragma once


#include "bertini2/system.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"

#include "bertini2/blackbox/config.hpp"


namespace bertini{
namespace blackbox{





template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> UserHomSpecifyStart(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return ZeroDimSpecifyTracker<start_system::User>(rt, ts...);
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeUserHom(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return UserHomSpecifyStart(rt, ts...);
}


} //ns blackbox
} //ns bertini
