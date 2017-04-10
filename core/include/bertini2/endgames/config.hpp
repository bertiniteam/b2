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
// Copyright(C) 2015, 2016, 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


/**
\file include/bertini2/endgames/config.hpp

\brief Configs and settings for endgames.
*/

#pragma once

namespace bertini{

	

	namespace tracking{

		namespace endgame{


		// some forward declarations
		template<typename Tracker, typename Enable = void>
		class FixedPrecPowerSeriesEndgame;

		template<typename Tracker, typename Enable = void>
		class FixedPrecCauchyEndgame;

		class AMPPowerSeriesEndgame;
		class AMPCauchyEndgame;

		template<typename TrackerType, typename FinalEGT> 
		class PowerSeriesEndgame;

		template<typename TrackerType, typename FinalEGT> 
		class CauchyEndgame;
		}

		
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
		{ };

		template<>
		struct EndgameSelector<DoublePrecisionTracker>
		{
			using PSEG = endgame::FixedPrecPowerSeriesEndgame<DoublePrecisionTracker>;
			using Cauchy = endgame::FixedPrecCauchyEndgame<DoublePrecisionTracker>;
		};

		template<>
		struct EndgameSelector<MultiplePrecisionTracker>
		{
			using PSEG = endgame::FixedPrecPowerSeriesEndgame<MultiplePrecisionTracker>;
			using Cauchy = endgame::FixedPrecCauchyEndgame<MultiplePrecisionTracker>;
		};

		template<class D>
		struct EndgameSelector<FixedPrecisionTracker<D> >
		{
			using PSEG = typename EndgameSelector<D>::PSEG;
			using Cauchy = typename EndgameSelector<D>::Cauchy;
		};

		template<>
		struct EndgameSelector<AMPTracker>
		{
			using PSEG = endgame::AMPPowerSeriesEndgame;
			using Cauchy = endgame::AMPCauchyEndgame;
		};

	namespace config{

		template<typename T>
		struct Security
		{
			int level = 0; //SecurityLevel
			T max_norm = T(100000); //SecurityMaxNorm wrong default value
		};

		template<typename T>
		struct Endgame
		{
			unsigned num_sample_points = 3; //NumSamplePoints default = 2
			T min_track_time = T(1e-100); //nbrh radius in Bertini book. NbhdRadius
			T sample_factor = T(1)/T(2); //SampleFactor
			unsigned max_num_newton_iterations = 15; // the maximum number allowable iterations during endgames, for points used to approximate the final solution.

			T final_tolerance = 1e-11;///< The tolerance to which to compute the endpoint using the endgame.
		};


		struct PowerSeries
		{
			unsigned max_cycle_number = 6; //MaxCycleNum
			unsigned cycle_number_amplification = 5;
		};

		template<typename T>
		struct Cauchy
		{
			T cycle_cutoff_time = T(1)/T(100000000); //CycleTimeCutoff
			T ratio_cutoff_time = T(1)/T(100000000000000); //RatioTimeCutoff
			T minimum_for_c_over_k_stabilization = T(3)/T(4);
			unsigned int num_needed_for_stabilization = 3;
			T maximum_cauchy_ratio = T(1)/T(2);
			unsigned int fail_safe_maximum_cycle_number = 250; //max number of loops before giving up. 

		};


		struct TrackBack
		{
			unsigned minimum_cycle; //MinCycleTrackback, default = 4
			bool junk_removal_test; //JunkRemovalTest, default = 1
			unsigned max_depth_LDT; //MaxLDTDepth, default = 3
		};

	}



	template <typename T>
	struct AlgoTraits;

	/**
	specialization for PowerSeries, which uses CRTP
	*/
	template<typename TrackerType, typename FinalEGT>
	struct AlgoTraits< typename endgame::PowerSeriesEndgame<TrackerType, FinalEGT>>
	{
		using NeededConfigs = detail::TypeList<
			config::PowerSeries,
			config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Security<typename TrackerTraits<TrackerType>::BaseRealType>
			>;
	};

	
	template<typename TrackerType>
	struct AlgoTraits<typename endgame::FixedPrecPowerSeriesEndgame<TrackerType>> : 
		public AlgoTraits<
				endgame::PowerSeriesEndgame<TrackerType, typename endgame::FixedPrecPowerSeriesEndgame<TrackerType>>
				>
	{
		using NeededConfigs = detail::TypeList<
			config::PowerSeries,
			config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Security<typename TrackerTraits<TrackerType>::BaseRealType>
			>;
	};

	template<>
	struct AlgoTraits<endgame::AMPPowerSeriesEndgame> : 
		public AlgoTraits<
				endgame::PowerSeriesEndgame<AMPTracker, typename endgame::AMPPowerSeriesEndgame>
				>
	{
		using NeededConfigs = detail::TypeList<
			config::PowerSeries,
			config::Endgame<typename TrackerTraits<AMPTracker>::BaseRealType>,
			config::Security<typename TrackerTraits<AMPTracker>::BaseRealType>
			>;
	};

	/**
	specialization for Cauchy, which uses CRTP
	*/
	template<typename TrackerType, typename FinalEGT>
	struct AlgoTraits< typename endgame::CauchyEndgame<TrackerType, FinalEGT>>
	{
		using NeededConfigs = detail::TypeList<
			config::Cauchy<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Security<typename TrackerTraits<TrackerType>::BaseRealType>>;
	};
	
	template<typename TrackerType>
	struct AlgoTraits<typename endgame::FixedPrecCauchyEndgame<TrackerType>> : 
		public AlgoTraits<
				endgame::CauchyEndgame<TrackerType, typename endgame::FixedPrecCauchyEndgame<TrackerType>>
				>
	{
		using NeededConfigs = detail::TypeList<
			config::Cauchy<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Endgame<typename TrackerTraits<TrackerType>::BaseRealType>,
			config::Security<typename TrackerTraits<TrackerType>::BaseRealType>>;
	};
	
	template<>
	struct AlgoTraits<endgame::AMPCauchyEndgame> : 
		public AlgoTraits<
				endgame::CauchyEndgame<AMPTracker, typename endgame::AMPCauchyEndgame>
				>
	{
		using NeededConfigs = detail::TypeList<
			config::Cauchy<typename TrackerTraits<AMPTracker>::BaseRealType>,
			config::Endgame<typename TrackerTraits<AMPTracker>::BaseRealType>,
			config::Security<typename TrackerTraits<AMPTracker>::BaseRealType>
			>;
	};



	}
}