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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Notre Dame
//  Summer 2016, Spring 2018
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
#include "generic_observable.hpp"
#include "configured_visitor.hpp"

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
			friend class ::boost::python::def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			using ComplexT = typename TrackerTraits<TrackerT>::BaseComplexT;
			using RealT = typename TrackerTraits<TrackerT>::BaseRealT;


			// resolve overloads for getting and setting predictor method.
			void (TrackerT::*set_predictor_)(Predictor)= &TrackerT::SetPredictor;
			Predictor (TrackerT::*get_predictor_)(void) const = &TrackerT::GetPredictor;
			
			// start_time and end_time are intentionally taken by value, not by const&.
			// eigenpy's from-python converter for the writable Eigen::Ref<Vec<ComplexT>>
			// 'result' argument corrupts the boost.python converter storage of an adjacent
			// scalar mpc_complex const& argument, leaving start_time/end_time as garbage
			// (invalid mpfr limb pointers -> MPFR set_prec assertion on first use).
			// Passing the scalars by value forces an independent copy that side-steps the
			// clobbered converter storage. See git history for the full diagnosis.
			static
			SuccessCode track_path_wrap(TrackerT const& self, Eigen::Ref<Vec<ComplexT>> result, ComplexT start_time, ComplexT end_time, Vec<ComplexT> const& start_point)
			{
				Vec<ComplexT> temp_result(self.GetSystem().NumVariables());
				auto code = self.TrackPath(temp_result, start_time, end_time, start_point);
				result = temp_result;
				return code;
			}



		};// TrackerVisitor class
		
		
		/**
		 AMP Tracker class
		 */
		template<typename TrackerT>
		class AMPTrackerVisitor: public def_visitor<AMPTrackerVisitor<TrackerT> >
		{
			friend class ::boost::python::def_visitor_access;
			
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
			friend class ::boost::python::def_visitor_access;
			
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
			friend class ::boost::python::def_visitor_access;
			
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
			friend class ::boost::python::def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		};// SteppingVisitor class



		
		// template<typename NumT>
		// class TolerancesVisitor: public def_visitor<TolerancesVisitor<NumT> >
		// {
		// 	friend class ::boost::python::def_visitor_access;

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
		
		
		// Visitor body definitions — template member functions in header so split TUs can instantiate.

		template<typename TrackerT>
		template<class PyClass>
		void TrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("setup", &TrackerT::Setup, (arg("predictor"), arg("tolerance"), arg("truncation"), arg("stepping"),arg("newton")), "Set values for the internal configuration of the tracker.  tolerance and truncation are both real doubles.  predictor is a valid value for predictor choice.  stepping and newton are the config structs from bertini.tracking.")

			.def("track_path", &track_path_wrap,
				 (arg("self"),arg("result"), "start_time", "end_time", "start_point"),
				 "The main function of the tracker, once its set up.  The first argument is the output.  Feed it, in (result, start_time, end_time, start_point")

			.def("get_system",&TrackerT::GetSystem,return_internal_reference<>(), "Gets an internal reference to the tracked system.")

			.def("predictor",get_predictor_,(arg("self")), "Query the current predictor method used by the tracker.")
			.def("predictor",set_predictor_,(arg("self"), arg("predictor")), "Set the predictor method used by the tracker.")

			.def("set_stepsize", &TrackerT::SetStepSize, (arg("self"), arg("stepsize")),"Set the stepsize for the tracker")

			.def("reinitialize_initial_step_size", &TrackerT::ReinitializeInitialStepSize, (arg("self"), arg("val")), "Set whether the tracker should re-set the stepsize to the configured-initial stepsize when it starts tracking.  Feed it a bool")
			.def("num_total_steps_taken", &TrackerT::NumTotalStepsTaken, (arg("self")),"Ask how many steps have been taken so far, including failures")

			.def("tracking_tolerance", &TrackerT::TrackingTolerance, (arg("self")), "Get.  A step is labeled as a failure if newton correcting doesn't yield a residual less than this tolerance.  A real number, the smaller the slower tracking, generally speaking")
			.def("tracking_tolerance", &TrackerT::SetTrackingTolerance, (arg("self"), arg("tol")), "Set the tracking tolerance for the tracker")

			.def("infinite_truncation_tolerance", &TrackerT::SetInfiniteTruncationTolerance, (arg("self"), arg("tol")) ,"Set the path truncation tolerance for infinite paths for the tracker")
			.def("infinite_truncation_tolerance", &TrackerT::InfiniteTruncationTolerance, (arg("self")), "Get the path truncation tolerance for infinite paths for the tracker")

			.def("infinite_truncation", &TrackerT::SetInfiniteTruncation, (arg("self"), arg("val")), "Decide whether the tracker should truncate infinite paths.  See also infinite_truncation_tolerance")
			.def("infinite_truncation", &TrackerT::InfiniteTruncation, (arg("self")), "Get the bool for whether the tracker should truncate infinite paths.  See also infinite_truncation_tolerance")

			.def("get_stepping",&TrackerT::template Get<tracking::SteppingConfig>,return_internal_reference<>(), (arg("self")), "Get the tracker's internal configuration for things that control stepping behaviour")
			.def("get_newton",&TrackerT::template Get<tracking::NewtonConfig>,return_internal_reference<>(), (arg("self")), "Get the tracker's internal configuration for Newton correction")
			.def("set_stepping",&TrackerT::template Set<tracking::SteppingConfig>, (arg("self"), arg("config")), "Set the tracker's internal configuration for things that control stepping behaviour")
			.def("set_newton",&TrackerT::template Set<tracking::NewtonConfig>, (arg("self"), arg("config")), "Set the tracker's internal configuration for Newton correction")

			.def("current_point", &TrackerT::CurrentPoint, (arg("self")), "what is the current point?")
			.def("current_time", &TrackerT::CurrentTime, (arg("self")), "what is the current time?")
			.def("current_precision", &TrackerT::CurrentPrecision, (arg("self")), "what is the current working precision?")

			// generic, type-list-driven config interface (get_config/set_config/config_types).
			// the get_stepping/set_stepping/get_newton/set_newton above remain as convenience aliases.
			.def(ConfiguredVisitor<TrackerT>())

			.def(ObservableVisitor<TrackerT>());
			;
		}


		template<typename TrackerT>
		template<class PyClass>
		void AMPTrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("precision_setup", &TrackerT::PrecisionSetup)
			.def("precision_preservation", &TrackerT::PrecisionPreservation, "Turn on or off the preservation of precision.  That is, if this is on (true), then the precision of the final point will be the precision of the start point.  Generally, you want to let precision drift, methinks.")

			.def("refine", return_Refine3_ptr<dbl>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			.def("refine", return_Refine3_ptr<mpfr_complex>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			.def("refine", return_Refine4_ptr<dbl>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time"), arg("tolerance"), arg("max_iterations")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			.def("refine", return_Refine4_ptr<mpfr_complex>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time"), arg("tolerance"), arg("max_iterations")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			;
		}


		template<typename TrackerT>
		template<class PyClass>
		void FixedDoubleTrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("refine", return_Refine3_ptr<dbl>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")

			.def("refine", return_Refine4_ptr<dbl>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time"), arg("tolerance"), arg("max_iterations")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			;
		}


		template<typename TrackerT>
		template<class PyClass>
		void FixedMultipleTrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("refine", return_Refine3_ptr<mpfr_complex>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")

			.def("refine", return_Refine4_ptr<mpfr_complex>(),
				(arg("self"), arg("result"), arg("start_point"), arg("time"), arg("tolerance"), arg("max_iterations")),
				"refine a point using this tracker, from `start_point`, at fixed `time`.  returns a success code, computed refined point is in `result`.")
			;
		}


		template<typename T>
		template<class PyClass>
		void SteppingVisitor<T>::visit(PyClass& cl) const
		{
			// initial_step_size, max_step_size, step_size_success_factor, step_size_fail_factor are
			// stored as mpq_rational (no MPFR precision state) but exposed to Python as mpfr_float
			// so the existing string/Float setter API is unchanged.  The round-trip is exact:
			// every mpfr_float has an exact rational representation, and converting back recovers it.
			cl
			.add_property("initial_step_size",
				+[](tracking::SteppingConfig const& c) -> mpfr_float { return mpfr_float(c.initial_step_size); },
				+[](tracking::SteppingConfig& c, mpfr_float const& v) { c.initial_step_size = mpq_rational(v); },
				"The initial stepsize when tracking is started.  See also tracking.AMPTracker.reinitialize_initial_step_size")
			.add_property("max_step_size",
				+[](tracking::SteppingConfig const& c) -> mpfr_float { return mpfr_float(c.max_step_size); },
				+[](tracking::SteppingConfig& c, mpfr_float const& v) { c.max_step_size = mpq_rational(v); },
				"The maximum allowed stepsize during tracking.  See also min_num_steps")
			.def_readwrite("min_step_size", &tracking::SteppingConfig::min_step_size,"The minimum stepsize the tracker is allowed to take.  See also max_step_size")
			.add_property("step_size_success_factor",
				+[](tracking::SteppingConfig const& c) -> mpfr_float { return mpfr_float(c.step_size_success_factor); },
				+[](tracking::SteppingConfig& c, mpfr_float const& v) { c.step_size_success_factor = mpq_rational(v); },
				"The scale factor for stepsize, after some consecutive steps.  See also consecutive_successful_steps_before_stepsize_increase")
			.add_property("step_size_fail_factor",
				+[](tracking::SteppingConfig const& c) -> mpfr_float { return mpfr_float(c.step_size_fail_factor); },
				+[](tracking::SteppingConfig& c, mpfr_float const& v) { c.step_size_fail_factor = mpq_rational(v); },
				"The scale factor for stepsize, after a fail happens.  See also step_size_success_factor")
			.def_readwrite("consecutive_successful_steps_before_stepsize_increase", &tracking::SteppingConfig::consecutive_successful_steps_before_stepsize_increase,"This number of successful steps have to taken consecutively, and then the stepsize is permitted to increase")
			.def_readwrite("min_num_steps", &tracking::SteppingConfig::min_num_steps, "The minimum number of steps the tracker can take between now and then.  This is useful if you are tracking closely between times, and want to guarantee some number of steps are taken.  Then again, this could be wasteful, too.")
			.def_readwrite("max_num_steps", &tracking::SteppingConfig::max_num_steps, "The maximum number of steps.  Tracking will die if it tries to take more than this number, sad day.")
			.def_readwrite("frequency_of_CN_estimation", &tracking::SteppingConfig::frequency_of_CN_estimation, "How frequently the condition number should be updated.  Less frequently is faster (estimation requires an additional linear solve), but may cause precision adjustment to lag behind.")
			;
		}


		// Prototypes for expose functions defined in split .cpp files.
		void ExportTrackers();
		void ExportAMPTracker();
		void ExportFixedTrackers();
		void ExportFixedDoubleTracker();
		void ExportFixedMultipleTracker();

		void ExportConfigSettings();

}}// re: namespaces




