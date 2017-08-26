//This file is part of Bertini 2.
//
//python/endgame_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Dani Brake
//  University of Notre Dame
//  Summer 2016
//
//
//  python/endgame_export.hpp:  Header file for exposing endgames to python.

#pragma once

#include "python_common.hpp"

#include <bertini2/endgames.hpp>

namespace bertini{
	namespace python{

		using namespace bertini::tracking;

		/**
		 Abstract Endgame class
		 */
		template<typename EndgameT>
		class EndgameVisitor: public def_visitor<EndgameVisitor<EndgameT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			using CT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseComplexType;
			using RT = typename TrackerTraits<typename EndgameT::TrackerType>::BaseRealType;

			template <typename T>
			using success_time_space = SuccessCode (EndgameT::*)(T const&, Vec<T> const&);
			template <typename T>
			using success_time_space_time = SuccessCode (EndgameT::*)(T const&, Vec<T> const&, T const&);

			template <typename T>
			static success_time_space<T> RunDefaultTime()
			{
				return &EndgameT::Run;
			};

			template <typename T>
			static success_time_space_time<T> RunCustomTime()
			{
				return &EndgameT::Run;
			};
			

			unsigned (EndgameT::*get_cycle_number_)() const = &EndgameT::CycleNumber;
		};// EndgameVisitor class



		/**
		 Particulars for the PowerSeries endgame.
		 */
		template<typename PowerSeriesT>
		class PowerSeriesVisitor: public def_visitor<PowerSeriesVisitor<PowerSeriesT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			using CT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseComplexType;
			using RT = typename TrackerTraits<typename PowerSeriesT::TrackerType>::BaseRealType;


		};// CauchyVisitor class




		/**
		 Particulars for the Cauchy endgame.
		 */
		template<typename CauchyT>
		class CauchyVisitor: public def_visitor<CauchyVisitor<CauchyT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			using CT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseComplexType;
			using RT = typename TrackerTraits<typename CauchyT::TrackerType>::BaseRealType;


		};// CauchyVisitor class





		class SecurityVisitor: public def_visitor<SecurityVisitor>
		{
			friend class def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def_readwrite("level", &endgame::SecurityConfig::level,"Turns on or off truncation of paths going to infinity during the endgame.  0 is off, 1 is on.")
				.def_readwrite("max_norm", &endgame::SecurityConfig::max_norm,"If on, the norm at which to truncate a path.")
				;
			}

		};




		
		class EndgameConfigVisitor: public def_visitor<EndgameConfigVisitor>
		{
			friend class def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def_readwrite("num_sample_points", &endgame::EndgameConfig::num_sample_points,"The number of points to use for extrapolant calculation.  In the Power Series Endgame, the is the number of geometrically spaces points on the path.  For Cauchy, this is the number of points on each circle tracked around the target time value.")
				.def_readwrite("min_track_time", &endgame::EndgameConfig::min_track_time,"The minimum distance from the target time to track to.  Decreasing this may help failing runs succeed, or maybe not, because you are, after all, tracking toward a singularity.")
				.def_readwrite("sample_factor", &endgame::EndgameConfig::sample_factor,"The factor by which to space the geometrically spaced `distance' between sample points, or sample circles for Cauchy.")
				.def_readwrite("max_num_newton_iterations", &endgame::EndgameConfig::max_num_newton_iterations,"the maximum number of newton iterations to be taken during sample point sharpening.  Increasing this can help speed convergence, at the risk of path jumping.")
				;
			}

		};


		class PowerSeriesConfigVisitor: public def_visitor<PowerSeriesConfigVisitor>
		{
			friend class def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def_readwrite("max_cycle_number", &endgame::PowerSeriesConfig::max_cycle_number,"The maximum cycle number to consider, when calculating the cycle number which best fits the path being tracked.")
				.def_readwrite("cycle_number_amplification", &endgame::PowerSeriesConfig::cycle_number_amplification,"The maximum number allowable iterations during endgames, for points used to approximate the final solution.")
				;
			}

		};





		class CauchyConfigVisitor: public def_visitor<CauchyConfigVisitor>
		{
			friend class def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def_readwrite("cycle_cutoff_time", &endgame::CauchyConfig::cycle_cutoff_time)
				.def_readwrite("ratio_cutoff_time", &endgame::CauchyConfig::ratio_cutoff_time)
				.def_readwrite("minimum_for_c_over_k_stabilization", &endgame::CauchyConfig::minimum_for_c_over_k_stabilization)
				.def_readwrite("maximum_cauchy_ratio", &endgame::CauchyConfig::maximum_cauchy_ratio)
				.def_readwrite("fail_safe_maximum_cycle_number", &endgame::CauchyConfig::fail_safe_maximum_cycle_number, "max number of loops before giving up." )
				;
			}

		};




		// now prototypes for expose functions defined in the .cpp files for the python bindings.

		/**
		The main function for exporting the bound endgames to Python.

		This should be the only function called from the main function defining the module, and should call all those functions exposing particular endgames.
		*/
		void ExportEndgames();


		/**
		export the power series endgame incarnations
		*/
		void ExportAMPPSEG();
		void ExportFDPSEG();
		void ExportFMPSEG();


		/**
		export the cauchy endgame incarnations
		*/
		void ExportAMPCauchyEG();
		void ExportFDCauchyEG();
		void ExportFMCauchyEG();

		

}}// re: namespaces

