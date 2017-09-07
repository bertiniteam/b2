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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


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
	\brief Converged -- endgame is done!
	*/
	ADD_BERTINI_EVENT_TYPE(Converged,EndgameSuccess);

	/**
	\brief Refined a sample
	*/
	ADD_BERTINI_EVENT_TYPE(SampleRefined,EndgameEvent);

	}// re: namespace endgames
}// re: namespace bertini

