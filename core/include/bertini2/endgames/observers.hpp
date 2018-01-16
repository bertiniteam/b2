//This file is part of Bertini 2.
//
//endgames/observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//endgames/observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with endgames/observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


/**
\file endgames/observers.hpp

\brief Contains the endgames/observers base types
*/

#pragma once

#include "bertini2/endgames/events.hpp"
#include "bertini2/logging.hpp"
#include "bertini2/detail/observer.hpp"

#include <boost/type_index.hpp>

namespace bertini {

	namespace endgame{


/**
\brief Logs the endgame run, with gory detail.

\ingroup loggers observers
*/
template <typename EndgameT>
struct GoryDetailLogger : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;
using BCT = typename EndgameT::BaseComplexType;

virtual ~GoryDetailLogger() = default;

virtual void Observe(AnyEvent const& e) override
{
	if(auto p = dynamic_cast<const TimeAdvanced<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "time advanced " << p->Get().LatestTime();
	}
	
	else if (auto p = dynamic_cast<const SampleRefined<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "refined a sample, huzzah";
	}

	else if (auto p = dynamic_cast<const CircleAdvanced<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "advanced around the circle, to " << p->NewSample()<< " at time " << p->NewTime();
	}

	else if (auto p = dynamic_cast<const ClosedLoop<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "closed a loop, cycle number " << p->Get().CycleNumber();
	}
	else if (auto p = dynamic_cast<const ApproximatedRoot<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "approximated the target root.  approximation " << p->Get().template FinalApproximation<BCT>() << " with error " << p->Get().ApproximateError();
	}

	else if (auto p = dynamic_cast<const PrecisionChanged<AMPEndgame>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "precision changed from  " << p->Previous() << " to " << p->Next();
	}

	else if (auto p = dynamic_cast<const InEGOperatingZone<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "made it to the endgame operating zone at time " << p->Get().LatestTime();
	}

	else if(auto p = dynamic_cast<const Converged<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "converged at time " << p->Get().LatestTime() << " with result " << p->Get().template FinalApproximation<BCT>() << " and residual " << p->Get().ApproximateError();
	}
	else if (auto p = dynamic_cast<const Initializing<EmitterT>*>(&e))
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "starting running " << boost::typeindex::type_id<EmitterT>().pretty_name();
	}
	else
	{
		BOOST_LOG_TRIVIAL(severity_level::debug) << "unprogrammed response for event of type " << boost::typeindex::type_id_runtime(e).pretty_name();
	}
}

}; // gory detail


	} //re: namespace endgames
}// re: namespace bertini
