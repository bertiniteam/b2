//This file is part of Bertini 2.
//
//bertini2/common/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/common/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/include/bertini2/trackers/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


#ifndef BERTINI2_COMMON_CONFIG
#define BERTINI2_COMMON_CONFIG

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
        /// No tracking was attempted.  This is a DEFAULT, never a verdict: no tracker returns it,
        /// and none may start to, because callers read it as "we never touched this path" and act
        /// on that -- a solve records nothing for such a path, counts it as unreached, and reports
        /// itself cut short.  A path that was attempted and got nowhere has a code that says what
        /// stopped it (SingularStartPoint from the initial refinement, MaxNumStepsTaken and its
        /// siblings from the per-iteration budgets, ExternallyTerminated, WallClockLimitReached).
        /// Pinned by test: a_tracker_never_returns_never_started.
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
        FailedToSelectPrecisionAndStepsize,
        WallClockLimitReached,   ///< The tracker's wall-clock deadline passed between steps; the path was abandoned where it was.  Append new values AFTER this one: the integers are part of the b2rec record contract.

    };

    // NumErrorT (the error/tolerance numeric type) now lives in num_traits.hpp -- the foundational
    // header included nearly everywhere -- so tolerance-typed parameters can name it without any
    // header adding a new include.  It is NOT redefined here to keep a single source of truth.
} // namespace bertini

#include "bertini2/common/stream_enum.hpp"


#endif // include guard
