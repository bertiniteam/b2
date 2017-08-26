//This file is part of Bertini 2.
//
//python/tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  James Collins
//  West Texas A&M University
//  Summer 2016
//
//
//
//  python/tracker.hpp:  Header file for exposing trackers to python.

#pragma once

#include "python_common.hpp"

#include <bertini2/trackers/tracker.hpp>

namespace bertini{
	namespace python{


		using namespace bertini::tracking;

		/**
		 Abstract Tracker class
		 */
		template<typename TrackerT>
		class TrackerVisitor: public def_visitor<TrackerVisitor<TrackerT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			using CT = typename TrackerTraits<TrackerT>::BaseComplexType;
			using RT = typename TrackerTraits<TrackerT>::BaseRealType;


			// resolve overloads for getting and setting predictor method.
			void (TrackerT::*set_predictor_)(Predictor)= &TrackerT::SetPredictor;
			Predictor (TrackerT::*get_predictor_)(void) const = &TrackerT::GetPredictor;
			





		};// TrackerVisitor class
		
		
		/**
		 AMP Tracker class
		 */
		template<typename TrackerT>
		class AMPTrackerVisitor: public def_visitor<AMPTrackerVisitor<TrackerT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			// resolve overloads for refining a point.
			template <typename T>
			using Refine3_ptr = SuccessCode (TrackerT::*)(Vec<T>&, Vec<T> const&, T const&) const;
			template <typename T>
			static Refine3_ptr<T> return_Refine3_ptr()
			{
				return &TrackerT::template Refine<T>;
			};
			
			template <typename ComplexT>
			using Refine4_ptr = SuccessCode (TrackerT::*)(Vec<ComplexT>&, Vec<ComplexT> const&, ComplexT const&, double const&, unsigned) const;
			template <typename ComplexT>
			static Refine4_ptr<ComplexT> return_Refine4_ptr()
			{
				return &TrackerT::template Refine<ComplexT>;
			};

			
		};// AMPTrackerVisitor class

		
		/**
		  Fixed Double Tracker class
		 */
		template<typename TrackerT>
		class FixedDoubleTrackerVisitor: public def_visitor<FixedDoubleTrackerVisitor<TrackerT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			// resolve overloads for refining a point.
			template <typename T>
			using Refine3_ptr = SuccessCode (TrackerT::*)(Vec<T>&, Vec<T> const&, T const&) const;
			template <typename T>
			static Refine3_ptr<T> return_Refine3_ptr()
			{
				return &TrackerT::template Refine<T>;
			};
			
			template <typename ComplexT>
			using Refine4_ptr = SuccessCode (TrackerT::*)(Vec<ComplexT>&, Vec<ComplexT> const&, ComplexT const&, double const&, unsigned) const;
			template <typename ComplexT>
			static Refine4_ptr<ComplexT> return_Refine4_ptr()
			{
				return &TrackerT::template Refine<ComplexT>;
			};

			
		};// FixedDoubleTrackerVisitor class

		
		/**
		 Fixed Multiple Tracker class
		 */
		template<typename TrackerT>
		class FixedMultipleTrackerVisitor: public def_visitor<FixedMultipleTrackerVisitor<TrackerT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			// resolve overloads for refining a point.
			template <typename T>
			using Refine3_ptr = SuccessCode (TrackerT::*)(Vec<T>&, Vec<T> const&, T const&) const;
			template <typename T>
			static Refine3_ptr<T> return_Refine3_ptr()
			{
				return &TrackerT::template Refine<T>;
			};
			
			template <typename ComplexT>
			using Refine4_ptr = SuccessCode (TrackerT::*)(Vec<ComplexT>&, Vec<ComplexT> const&, ComplexT const&, double const&, unsigned) const;
			template <typename ComplexT>
			static Refine4_ptr<ComplexT> return_Refine4_ptr()
			{
				return &TrackerT::template Refine<ComplexT>;
			};
			
			
		};// FixedMultipleTrackerVisitor class

		

		
		
		/**
		 Stepping struct
		 */
		template<typename T>
		class SteppingVisitor: public def_visitor<SteppingVisitor<T> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};// SteppingVisitor class



		
		// template<typename NumT>
		// class TolerancesVisitor: public def_visitor<TolerancesVisitor<NumT> >
		// {
		// 	friend class def_visitor_access;

		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const
		// 	{
		// 		cl
		// 		.def_readwrite("newton_before_endgame", &Tolerances<NumT>::newton_before_endgame)
		// 		.def_readwrite("newton_during_endgame", &Tolerances<NumT>::newton_during_endgame)
		// 		.def_readwrite("final_tolerance", &Tolerances<NumT>::final_tolerance)
		// 		.def_readwrite("final_tolerance_multiplier", &Tolerances<NumT>::final_tolerance_multiplier)
		// 		.def_readwrite("path_truncation_threshold", &Tolerances<NumT>::path_truncation_threshold)
		// 		.def_readwrite("final_tolerance_times_final_tolerance_multiplier", &Tolerances<NumT>::final_tolerance_times_final_tolerance_multiplier)
		// 		;
		// 	}

		// };
		
		
		// now prototypes for expose functions defined in the .cpp files for the python bindings.
		void ExportTrackers();
		void ExportAMPTracker();
		void ExportFixedTrackers();
		void ExportFixedDoubleTracker();
		void ExportFixedMultipleTracker();
		
		void ExportConfigSettings();

}}// re: namespaces




