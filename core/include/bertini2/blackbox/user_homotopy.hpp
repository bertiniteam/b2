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
// Copyright(C) Bertini2 Development Team
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
#include "bertini2/blackbox/switches_zerodim.hpp"   // ZeroDimRT


namespace bertini{
namespace blackbox{


// User homotopies use the RefToGiven policy (the user owns target/start/homotopy), so they
// get their own tracker/endgame dispatch chain -- distinct from the generic ZeroDimSpecify*
// ladder, which always clones.  ts... = (target, start_system, homotopy).

template <typename TrackerType, typename EndgameType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> UserHomSpecifyComplete(ConstTs const& ...ts)
{
	return std::make_unique<
			algorithm::ZeroDim<TrackerType, EndgameType, System, policy::RefToGiven>
			>(ts...);
}

template <typename TrackerType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> UserHomSpecifyEndgame(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.endgame)
	{
		case type::Endgame::PowerSeries:
			return UserHomSpecifyComplete<TrackerType, typename endgame::EndgameSelector<TrackerType>::PSEG>(ts...);
		case type::Endgame::Cauchy:
			return UserHomSpecifyComplete<TrackerType, typename endgame::EndgameSelector<TrackerType>::Cauchy>(ts...);
	}
	throw std::runtime_error("unrecognized endgame type in UserHomSpecifyEndgame");
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> UserHomSpecifyTracker(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.tracker)
	{
		case type::Tracker::FixedDouble:
			return UserHomSpecifyEndgame<tracking::DoublePrecisionTracker>(rt, ts...);
		case type::Tracker::FixedMultiple:
			return UserHomSpecifyEndgame<tracking::MultiplePrecisionTracker>(rt, ts...);
		case type::Tracker::Adaptive:
			return UserHomSpecifyEndgame<tracking::AMPTracker>(rt, ts...);
	}
	throw std::runtime_error("unrecognized tracker type in UserHomSpecifyTracker");
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeUserHom(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return UserHomSpecifyTracker(rt, ts...);
}


} //ns blackbox
} //ns bertini
