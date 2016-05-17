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
//  Daniel Brake
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

#include <bertini2/tracking/tracker.hpp>

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
			void (TrackerT::*set_predictor_)(config::Predictor)= &TrackerT::Predictor;
			config::Predictor (TrackerT::*get_predictor_)(void) const = &TrackerT::Predictor;
			
			// resolve overloads for refining a point.
			template <typename T>
			using Refine3_ptr = SuccessCode (TrackerT::*)(Vec<T>&, Vec<T> const&, T const&) const;
			template <typename T>
			static Refine3_ptr<T> return_Refine3_ptr()
			{
				return &TrackerT::template Refine<T>;
			};
			
			template <typename ComplexT, typename RealT>
			using Refine4_ptr = SuccessCode (TrackerT::*)(Vec<ComplexT>&, Vec<ComplexT> const&, ComplexT const&, RealT const&, unsigned) const;
			template <typename ComplexT, typename RealT>
			static Refine4_ptr<ComplexT, RealT> return_Refine4_ptr()
			{
				return &TrackerT::template Refine<ComplexT, RealT>;
			};





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
			
		};// AMPTrackerVisitor class



		
//		/**
//		 TrackerTrait struct
//		 */
//		template<typename TrackerT>
//		class TrackerTraitVisitor: public def_visitor<TrackerTraitVisitor<TrackerT> >
//		{
//			friend class def_visitor_access;
//			
//		public:
//			template<class PyClass>
//			void visit(PyClass& cl) const;
//			
//			
//			
//		};// TrackerTraitVisitor class

		
		
		// now prototypes for expose functions defined in the .cpp files for the python bindings.
		void ExportTrackers();
		void ExportAMPTracker();
		void ExportFixedTrackers();
		void ExportFixedDoubleTracker();
		void ExportFixedMultipleTracker();
		
		void ExportConfigSettings();

}}// re: namespaces




