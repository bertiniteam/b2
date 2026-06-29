//This file is part of Bertini 2.
//
//nag_algorithms/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nag_algorithms/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nag_algorithms/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Tim Hodges, Colorado State University
// jeb collins, west texas a&m
 
#pragma once

#include "bertini2/system/start_systems.hpp"

#include <type_traits>

namespace bertini{
	namespace algorithm{

template<typename T>
using SolnCont = std::vector<T>;

namespace classic{

enum class EndgameChoice
{
	PowerSeries = 1,
	Cauchy = 2
};

struct EndgameChoiceConfig
{
	EndgameChoice endgame = EndgameChoice::Cauchy;
};

enum class AlgoChoice
{
	EvalFunctions = -4,
	EvalFunctionJacobian = -3,
	NewtonIteration = -2,
	NewtonIterationCondNum = -1,
	ZeroDim = 0,
	NID = 1,
	SampleComponent = 2,
	MembershipTest = 3,
	ExtractWitnessSet = 4,
	WitnessSetProjection = 5,
	IsosingularStab = 6
};

} // namespace classic



struct TolerancesConfig
{	
	using T = NumErrorT;

	T newton_before_endgame = T(1)/T(100000); //E.4.1
	T newton_during_endgame = T(1)/T(1000000); //E.4.2

	T final_tolerance = T(1)/T(100000000000); //E.5.1

	T path_truncation_threshold = T(100000); //E.4.13
};
		
	
struct MidPathConfig
{
	using T = NumErrorT;

	T same_point_tolerance = T(1)/T(100000);
};



struct AutoRetrackConfig
{
	using T = NumErrorT;

	T midpath_decrease_tolerance_factor = T(1)/T(2);
};



struct SharpeningConfig
{
	using T = NumErrorT;

	unsigned sharpendigits; ///< how many digits should be correct after sharpening.
	
	// std::function<Vec<T>> sharpen_method_; ///< function taking a vector, and sharpening it.

	T function_residual_tolerance = Eigen::NumTraits<T>::dummy_precision(); ///< A polynomial (or any function, really) evaluated at a point is considered to be 0 if the magnitude is smaller than this value.  See also RatioTolerance.  **Note that this value depends on the current default precision when this scruct is constructed.**

	T ratio_tolerance = T(99)/T(100); ///<  A computed value is considered to be zero if the ratio of two different approximations is smaller than this value.  See also FunctionTolerance
};



struct RegenerationConfig
{
	using T = NumErrorT;

	bool remove_infinite_endpoints = true; ///<  Bool indicating whether endpoints during the regeneration start point buildup step which are infinite should be discarded.  If you are not interested in infinite solutions, ensure this is true.  RegenRemoveInf

	bool higher_dimension_check = true; ///< RegenHigherDimTest
	unsigned start_level = 0;
	// These are the SLICE-moving tracking tolerances (Bertini 1's SliceTol* family) -- the tolerances
	// for moving the linear slices during regeneration, kept separate from the main tracking
	// tolerances in TolerancesConfig.  The slice_ prefix makes every config field name unique across
	// structs, which is what lets a field be set on an owner without naming its struct
	// (owner.update(field=...) routes by field).
	T slice_newton_before_endgame; ///< Slice-moving tracking tolerance before the endgame.  SliceTolBeforeEG
	T slice_newton_during_endgame; ///< Slice-moving tracking tolerance during the endgame.  SliceTolDuringEG
	T slice_final_tolerance; ///< Final tolerance to track the slice move to, using the endgame.  SliceFinalTol
};



struct PostProcessingConfig{
	using T = NumErrorT;
	
	T real_threshold = T(1)/T(100000000); ///< Bertini 1's `ImagThreshold`.  Threshold on the imaginary part of a (dehomogenized) solution coordinate being 0: a point is real if the infinity norm of the imaginary parts is below this.  If the imag part exceeds this, the point is considered complex.  Currently, this is the implemented available way in Bertini2 for determining this, but there are other methods.  Smale's alpha theory provides ways to prove that a point is real.  If this is something you need, please consider adding the method to the library, for all to use!  Or, if this is technically beyond your C++ capabilities, add as an issue on the github page, and indicate it as a feature request.  B1 default 1e-8.

	T endpoint_finite_threshold = T(100000);  ///< Bertini 1's `EndpointFiniteThreshold`.  An endpoint is considered to be at infinity if the infinity norm of its *dehomogenized* coordinates is larger than this value.  This uses the same dehomogenize-then-infinity-norm computation the endgame uses for its `Security::max_norm` divergence check (a separate, smaller threshold for bailing out *during* the endgame).  There is also `path_truncation_threshold` in Tolerances, which tells the path tracker to die if exceeded.  B1 default 1e5.

	T same_point_tolerance_multiplier {T(10)}; ///< Bertini 1's `EndpointSameThreshold`.  A multiplier (>= 1) on `final_tolerance`: two endpoints are considered the same point if the infinity norm of the difference of their *dehomogenized* coordinates is below `final_tolerance * same_point_tolerance_multiplier`.  Keeping it a multiplier (rather than an absolute tolerance) means the same-point test always stays a fixed factor looser than the accuracy you tracked to, even if `final_tolerance` is changed.  B1 default 10.

	T condition_number_threshold {T(100000000)}; ///< Bertini 1's `CondNumThreshold`.  An endpoint is considered singular if it is the endpoint of multiple paths (multiplicity > 1), or if the approximation of the condition number (in the spectral norm, as estimated by the tracker) is larger than this value.  B1 default 1e8.
};

/**
\brief The default ambient precision for a ZeroDimConfig, by complex type.

Double-precision solves work at double; multiprecision solves work at the current default precision
-- which is the precision the MultiplePrecisionTracker is built at, so that in a fixed-multiple solve
the start points, the tracker, and the working ("thread") precision are all one and the same value.
(Defaulting to DoublePrecision() here left fixed-multiple solves with the tracker at DefaultPrecision()
but the ambient/thread precision at double -- a mismatch the tracker rejects at start.)
*/
template<typename ComplexT>
inline unsigned DefaultInitialAmbientPrecision()
{
	if constexpr (std::is_same<ComplexT, complex_dbl>::value)
		return DoublePrecision();
	else
		return DefaultPrecision();
}

// Not templated on the complex type: the three times are stored as exact, precision-free
// mpq_rational (real) and converted to the tracking complex type at the current working precision at
// use -- the same pattern SteppingConfig uses for its step sizes.  This keeps the config precision-
// agnostic (one struct, not a DoublePrec/Multiprec pair), and is in fact more correct under adaptive
// precision: a value like the endgame boundary 1/10 is materialized at full working precision rather
// than frozen at whatever precision the config happened to be constructed at.
//
// The times are real because a zero-dim solve tracks along the real t-axis (1 -> 0).  Complex
// homotopy times remain available at the endgame level (set directly), just not through this config.
struct ZeroDimConfig
{
	// Per-complex-type default; the ZeroDim algorithm overwrites this in DefaultSettingsSetup with
	// DefaultInitialAmbientPrecision<BaseComplexT>() (it knows its tracking type, this struct does not).
	unsigned initial_ambient_precision = DefaultPrecision();
	unsigned max_num_crossed_path_resolve_attempts = 2; ///< The maximum number of times to attempt to re-solve crossed paths at the endgame boundary.

	/// Number of worker threads for a shared-memory (non-MPI) solve.  0 = auto
	/// (std::thread::hardware_concurrency); 1 = serial (no thread pool).  Overridden by the
	/// OMP_NUM_THREADS environment variable when set.  See parallel::EffectiveThreadCount.
	/// Under MPI the per-rank thread count comes from OMP_NUM_THREADS, not this field.
	unsigned num_threads = 0;

	mpq_rational start_time{1};          ///< Homotopy start time (t=1).
	mpq_rational endgame_boundary{1, 10}; ///< Time at which tracking hands off to the endgame (t=1/10).
	mpq_rational target_time{0};         ///< Homotopy target time (t=0).

	std::string path_variable_name = "ZERO_DIM_PATH_VARIABLE";
};

struct MetaConfig
{
	classic::AlgoChoice tracktype = classic::AlgoChoice::ZeroDim;
};

/**
Global RNG seed for reproducible runs.  random_seed == 0 (the default) draws from
std::random_device and reports the effective seed so the run can be reproduced.
Set via `randomseed: N;` in the classic Bertini input file or bertini.set_random_seed(N)
in Python.  Must be applied before system construction (gamma, patch, TD-constants).
*/
struct RandomConfig
{
	unsigned long random_seed = 0;
};


// a forward declare
template <typename T>
	struct AlgoTraits;


} } // namespaces
