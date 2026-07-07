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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


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

\ingroup observer
*/
template <typename EndgameT>
struct GoryDetailLogger : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;  ///< The event-emitter type.
using BCT = typename EndgameT::BaseComplexT;  ///< The complex number type.

virtual ~GoryDetailLogger() = default;

virtual ObserveResult Observe(AnyEvent const& e) override
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

	return ObserveResult::KeepObserving;
}

}; // gory detail



/**
\brief Counts endgame events and captures their payloads, for tests and diagnostics.

Records how many of each event type were delivered, and captures the most-recent CircleAdvanced
point/time and the Converged approximation -- reading each payload (and the state accessors LatestTime /
FinalApproximation / ApproximateError) the way a real observer would.  This deliberately exercises the
full event-delivery path, including the numeric-type conversion the adaptive-numeric-type endgame
performs when it emits from its hardware-complex_dbl fast lane, and the slot-aware state accessors.

\ingroup observer
*/
template <typename EndgameT>
struct EventRecorder : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;                     ///< The endgame type emitting the observed events.
using BCT = typename EndgameT::BaseComplexT;    ///< The boundary complex type of the observed endgame.

unsigned num_events            = 0;            ///< Total number of events observed.
unsigned num_time_advanced     = 0;            ///< Number of TimeAdvanced events observed.
unsigned num_sample_refined    = 0;            ///< Number of SampleRefined events observed.
unsigned num_circle_advanced   = 0;            ///< Number of circle-advanced events observed.
unsigned num_closed_loop       = 0;            ///< Number of closed-loop events observed.
unsigned num_approximated_root = 0;            ///< Number of approximated-root events observed.
unsigned num_in_eg_zone        = 0;            ///< Number of in-endgame-zone events observed.
unsigned num_converged         = 0;            ///< Number of converged events observed.
unsigned num_precision_changed = 0;            ///< Number of precision-changed events observed.

Vec<BCT> last_circle_point;                    ///< The most-recent circle sample point captured from an event.
BCT      last_circle_time;                     ///< The most-recent circle sample time captured from an event.
Vec<BCT> converged_point;                      ///< The converged root point captured from a converged event.

virtual ObserveResult Observe(AnyEvent const& e) override
{
	++num_events;
	if (auto p = dynamic_cast<const TimeAdvanced<EmitterT>*>(&e))
	{ ++num_time_advanced; (void)p->Get().LatestTime(); }

	else if (dynamic_cast<const SampleRefined<EmitterT>*>(&e))
	{ ++num_sample_refined; }

	else if (auto p = dynamic_cast<const CircleAdvanced<EmitterT>*>(&e))
	{ ++num_circle_advanced; last_circle_point = p->NewSample(); last_circle_time = p->NewTime(); }

	else if (dynamic_cast<const ClosedLoop<EmitterT>*>(&e))
	{ ++num_closed_loop; }

	else if (auto p = dynamic_cast<const ApproximatedRoot<EmitterT>*>(&e))
	{ ++num_approximated_root; (void)p->Get().template FinalApproximation<BCT>(); (void)p->Get().ApproximateError(); }

	else if (dynamic_cast<const InEGOperatingZone<EmitterT>*>(&e))
	{ ++num_in_eg_zone; }

	else if (auto p = dynamic_cast<const Converged<EmitterT>*>(&e))
	{ ++num_converged; converged_point = p->Get().template FinalApproximation<BCT>(); (void)p->Get().LatestTime(); }

	else if (dynamic_cast<const PrecisionChanged<AMPEndgame>*>(&e))
	{ ++num_precision_changed; }

	return ObserveResult::KeepObserving;
}

}; // EventRecorder


	} //re: namespace endgames
}// re: namespace bertini
