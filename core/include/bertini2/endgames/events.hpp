//This file is part of Bertini 2.
//
//events.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//events.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with events.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


/**
\file endgames/events.hpp

\brief Contains the endgames/events base types
*/

#pragma once

#include "bertini2/detail/events.hpp"

namespace bertini {

	namespace endgame{


	
	/**
	\brief Generic event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameEvent,ConstEvent);

	/**
	\brief Generic failure event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameFailure,ConstEvent);

	/**
	\brief Generic success event for Endgames
	*/
	ADD_BERTINI_EVENT_TYPE(EndgameSuccess,ConstEvent);

	/**
	\brief MinTrackTime reached during endgame
	*/
	ADD_BERTINI_EVENT_TYPE(MinTrackTimeReached,EndgameFailure);

	/**
	\brief Refining a sample failed for some reason
	*/
	ADD_BERTINI_EVENT_TYPE(RefiningFailed,EndgameFailure);

	/**
	\brief Cycle number computed was too high
	*/
	ADD_BERTINI_EVENT_TYPE(CycleNumTooHigh,EndgameFailure);

	/**
	\brief Security max norm reached. 
	*/
	ADD_BERTINI_EVENT_TYPE(SecurityMaxNormReached,EndgameFailure);

	/**
	\brief Started running the endgame
	*/
	ADD_BERTINI_EVENT_TYPE(Initializing,EndgameEvent);


	/**
	\brief Time advancing
	*/
	ADD_BERTINI_EVENT_TYPE(TimeAdvanced,EndgameEvent);

	/**
	\brief Advanced around the circle around target time
	*/
	template<class ObservedT>
	class CircleAdvanced : public EndgameEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:

		using CT = typename ObservedT::BaseComplexType;
		/**
		\brief The constructor for a CircleAdvanced Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		CircleAdvanced(const ObservedT & obs, 
		               Vec<CT> const& new_point,
		               CT const& new_time) : EndgameEvent<ObservedT>(obs),
													new_point_(new_point),
													new_time_(new_time)
		{}


		virtual ~CircleAdvanced() = default;
		CircleAdvanced() = delete;
		
		const auto& NewSample() const {return new_point_;}

		const auto& NewTime() const {return new_time_;}

	private:
		const Vec<CT>& new_point_;
		const CT& new_time_;
	};


	/**
	\brief Walked a complete loop around the target time.
	*/
	ADD_BERTINI_EVENT_TYPE(ClosedLoop,EndgameEvent);

	/**
	\brief Approximated a root at target time.
	*/
	ADD_BERTINI_EVENT_TYPE(ApproximatedRoot,EndgameEvent);

	/**
	\brief Converged -- endgame is done!
	*/
	ADD_BERTINI_EVENT_TYPE(Converged,EndgameSuccess);

	/**
	\brief Refined a sample
	*/
	ADD_BERTINI_EVENT_TYPE(SampleRefined,EndgameEvent);

	/**
	\brief Made it into the EG operating zone, or so we believe
	*/
	ADD_BERTINI_EVENT_TYPE(InEGOperatingZone,EndgameEvent);
	

	template<class ObservedT>
	class PrecisionChanged : public EndgameEvent<ObservedT>
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		/**
		\brief The constructor for a PrecisionChanged Event.

		\param obs The observable emitting the event.
		\param previous The precision before changing.
		\param next The precision after changing.
		*/
		PrecisionChanged(const ObservedT & obs, 
		             unsigned previous, unsigned next) : EndgameEvent<ObservedT>(obs),
													prev_(previous),
													next_(next)
		{}


		virtual ~PrecisionChanged() = default;
		PrecisionChanged() = delete;
		
		/**
		\brief Get the previous precision.
		*/
		auto Previous() const {return prev_;}

		/**
		\brief Get the next precision, what it changed to.
		*/
		auto Next() const {return next_;}
	private:
		const unsigned prev_, next_;
	};


	}// re: namespace endgames
}// re: namespace bertini

