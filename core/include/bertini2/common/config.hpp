//This file is part of Bertini 2.
//
//trackers/include/bertini2/trackers/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//trackers/include/bertini2/trackers/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/include/bertini2/trackers/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


// individual authors of this file include:
// dani brake, university of notre dame, university of wisconsin eau claire
// Tim Hodges, Colorado State University

#ifndef BERITNI2_COMMON_CONFIG
#define BERITNI2_COMMON_CONFIG

#pragma once


namespace bertini
{

	// aliases for the types used to contain space and time samples, and random vectors for the endgames.
	template<typename T> using SampCont = std::deque<Vec<T> >;
	template<typename T> using TimeCont = std::deque<T>;
	
	enum class ContStart{
		Front,
		Back
	};

	enum class SuccessCode
	{
		NeverStarted = -1,
		Success = 0,
		HigherPrecisionNecessary,
		ReduceStepSize,
		GoingToInfinity,
		FailedToConverge,
		MatrixSolveFailure,
		MatrixSolveFailureFirstPartOfPrediction,
		MaxNumStepsTaken,
		MaxPrecisionReached,
		MinStepSizeReached,
		Failure,
		SingularStartPoint,
		ExternallyTerminated,
		MinTrackTimeReached,
		SecurityMaxNormReached,
		CycleNumTooHigh,

	};

	using NumErrorT = double;
} // namespace bertini




#endif // include guard

