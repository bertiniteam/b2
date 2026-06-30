//This file is part of Bertini 2.
//
//include/bertini2/endgames/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/endgames/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/endgames/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


/**
\file include/bertini2/endgames/config.hpp

\brief Configs and settings for endgames.
*/

#pragma once

#include "bertini2/detail/typelist.hpp"
#include "bertini2/common/config.hpp"
#include "bertini2/mpfr_extensions.hpp"


namespace bertini{ namespace endgame{


		// // some forward declarations
// base
	template<class FlavorT, class PrecT>
	class EndgameBase;

// flavors
	template<typename PrecT> 
	class PowerSeriesEndgame;

	template<typename PrecT> 
	class CauchyEndgame;

// precision types
	template<typename TrackerT>
	class FixedPrecEndgame;

	class AMPEndgame;

	// end forward declarations



	/**
	Base empty class
	*/
	template <typename TrackerT>
	struct EGPrecSelector;

	// specialize this in the specific files it goes in, please

	



	/**
	\brief Facilitates lookup of required endgame type based on tracker type
	
	Your current choices are PSEG or Cauchy.

	To get the Power Series Endgame for Adaptive Precision Tracker, use the following example code:
	\code
	using EGT = EndgameSelector<AMPTracker>::PSEG
	\endcode

	\tparam TrackerT The type of tracker you want to use.
	*/
	template<typename TrackerT>
	struct EndgameSelector
	{ 
		using EGPrecT = typename EGPrecSelector<TrackerT>::type;  ///< The precision policy type for the tracker.

		using PSEG = endgame::PowerSeriesEndgame<EGPrecT>;   ///< The power-series endgame for this tracker.
		using Cauchy = endgame::CauchyEndgame<EGPrecT>;      ///< The Cauchy endgame for this tracker.
	};

	/// \brief Security (divergence-bailout) configuration for an endgame.
	struct SecurityConfig
	{
		int level = 0; ///< Security level; >0 disables the max-norm divergence check.
		NumErrorT max_norm = NumErrorT(1e4); ///< A path diverges if its norm exceeds this during the endgame.
	};


	/// \brief Configuration common to all endgames (sampling and final-approximation tolerances).
	struct EndgameConfig
	{
		using T = NumErrorT;  ///< The numeric (error) type.
		T sample_point_refinement_factor = 1e-2; ///< Extra tolerance for refining sample points before computing the final approximation.
		unsigned num_sample_points = 3; ///< Number of sample points used per approximation.
		T min_track_time = T(1e-100); ///< Smallest time the endgame will track to (the neighborhood radius).

		mpq_rational sample_factor{1, 2}; ///< Geometric factor between successive sample times (exact rational).

		unsigned max_num_newton_iterations = 15; ///< Maximum Newton iterations when refining endgame sample points.

		T final_tolerance = 1e-11;///< The tolerance to which to compute the endpoint using the endgame.

		// When the adaptive-numeric-type (double-first) endgame crosses a sample set from double up into
		// mpfr, the retained samples are widened (zero-padded) -- they were already tracked/refined to
		// final_tolerance, so by default we do NOT spend a Newton refine to sharpen them.  Set true to
		// refine every retained sample to the new (higher) precision immediately after a precision
		// increase.  No effect on fixed-precision endgames.
		bool refine_when_increasing_precision = false; ///< When true, re-refine retained samples after the endgame migrates to higher precision (default false).
	};


	/// \brief Power-series-endgame configuration.
	struct PowerSeriesConfig
	{
		unsigned max_cycle_number = 6; ///< Largest cycle number to consider.
		unsigned cycle_number_amplification = 5; ///< Multiplier bounding the search for the cycle number.
	};


	/// \brief Cauchy-endgame configuration.
	struct CauchyConfig
	{
		using T = NumErrorT;  ///< The numeric (error) type.

		T cycle_cutoff_time = T(1)/T(100000000); ///< Time below which the cycle-number heuristic stops.
		T ratio_cutoff_time = T(1)/T(100000000000000); ///< Time below which the c/k ratio test stops.
		T minimum_for_c_over_k_stabilization = T(3)/T(4); ///< Minimum c/k ratio accepted as stabilized.
		unsigned int num_needed_for_stabilization = 3; ///< Consecutive samples needed for c/k stabilization.
		T maximum_cauchy_ratio = T(1)/T(2); ///< Maximum accepted Cauchy ratio.
		unsigned int fail_safe_maximum_cycle_number = 250; ///< Max number of loops before giving up.

		// Number of consecutive circle-tracked Cauchy approximations that must report the SAME cycle
		// number before a converged approximation is trusted.  Guards against an UNRELIABLE cycle
		// number: when the working precision is too low to close the loop accurately (the circle of
		// radius |t| has enough path variation that tracking error breaks closure -- NOT monodromy:
		// at high precision a nonsingular endpoint closes at cycle 1 for every radius), the
		// loop-closing count thrashes (e.g. 41, 14, 36).  Requiring N consecutive identical cycle
		// numbers refuses to accept convergence until the estimate has genuinely settled.
		// See z_notes/20260629_endgame_stepsize_reset_rootcause.
		unsigned int num_consecutive_same_cycle_number = 2; ///< Consecutive identical cycle-number estimates required before trusting convergence (guards against an unreliable cycle number at too-low precision).

	};


	/// \brief Track-back endgame configuration.
	struct TrackBackConfig
	{
		unsigned minimum_cycle = 4; ///< Minimum cycle number for track-back.
		bool junk_removal_test = 1; ///< Whether to run the junk-removal test.
		unsigned max_depth_LDT = 3; ///< Maximum depth of the local dimension test.
	};


	// an empty base class
	template <typename T>
	struct AlgoTraits;

	/**
	specialization for PowerSeries, which uses CRTP
	*/
	template<typename PrecT>
	struct AlgoTraits< PowerSeriesEndgame<PrecT>>
	{
		/// The config types this endgame reads.
		using NeededConfigs = detail::TypeList<
			PowerSeriesConfig,
			EndgameConfig,
			SecurityConfig
			>;

		using EmitterType = PowerSeriesEndgame<PrecT>;  ///< The event-emitter type for this endgame.
	};


	/**
	specialization for Cauchy, which uses CRTP
	*/
	template<typename PrecT>
	struct AlgoTraits< CauchyEndgame<PrecT>>
	{
		/// The config types this endgame reads.
		using NeededConfigs = detail::TypeList<
			CauchyConfig,
			EndgameConfig,
			SecurityConfig>;
	};
	
	


} } // namespaces


