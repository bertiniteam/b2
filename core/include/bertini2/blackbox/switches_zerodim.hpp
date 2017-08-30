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
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/blackbox/switches_zerodim.hpp 

\brief A sequence of switches for getting a particular instantiation of a ZeroDim algorithm, based on runtime options.
*/



#pragma once


#include "bertini2/system.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"

namespace bertini{
namespace blackbox{

namespace type{
enum class Start{ TotalDegree, MHom, User};
enum class Tracker{ FixedDouble, FixedMultiple, Adaptive};
enum class Endgame{ PowerSeries, Cauchy};
}

struct ZeroDimRT
{
	type::Start start = type::Start::TotalDegree;
	type::Tracker tracker = type::Tracker::Adaptive;
	type::Endgame endgame = type::Endgame::Cauchy;
};


template <typename StartType, typename TrackerType, typename EndgameType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecify(ConstTs const& ...ts)
{
	return std::make_unique<
			algorithm::ZeroDim<
				TrackerType, 
				EndgameType, 
				System, 
				StartType>
			>(ts...);

}

template <typename StartType, typename TrackerType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyEndgame(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.endgame)
	{
		case type::Endgame::PowerSeries:
			return ZeroDimSpecify<StartType, TrackerType, 
					typename endgame::EndgameSelector<TrackerType>::PSEG>(ts...);

		case type::Endgame::Cauchy:
			return ZeroDimSpecify<StartType, TrackerType, 
					typename endgame::EndgameSelector<TrackerType>::Cauchy>(ts...);
	}
}

template <typename StartType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyTracker(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.tracker)
	{
		case type::Tracker::FixedDouble:
			return ZeroDimSpecifyEndgame<StartType, tracking::DoublePrecisionTracker>(rt, ts...);
		case type::Tracker::FixedMultiple:
			return ZeroDimSpecifyEndgame<StartType, tracking::MultiplePrecisionTracker>(rt, ts...);
		case type::Tracker::Adaptive:
			return ZeroDimSpecifyEndgame<StartType, tracking::AMPTracker>(rt, ts...);
	}
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyStart(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.start)
	{
		case type::Start::TotalDegree:
			return ZeroDimSpecifyTracker<start_system::TotalDegree>(rt, ts...);
		case type::Start::MHom:
			return ZeroDimSpecifyTracker<start_system::MHomogeneous>(rt, ts...);
		case type::Start::User:
			return ZeroDimSpecifyTracker<start_system::User>(rt, ts...);
	}
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeZeroDim(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return ZeroDimSpecifyStart(rt, ts...);
}


} //ns blackbox
} //ns bertini
