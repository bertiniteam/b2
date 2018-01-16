//This file is part of Bertini 2.
//
//python/tracker_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/tracker_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/tracker_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Danielle Brake
//  University of Notre Dame
//  Summer 2016, Spring 2018
//
//
//  python/tracker_export.cpp:  source file for exposing trackers to python.

#include "tracker_export.hpp"

namespace bertini{
	namespace python{

		template<typename TrackerT>
		template<class PyClass>
		void TrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("setup", &TrackerT::Setup, (boost::python::arg("predictor"), boost::python::arg("tolerance"), boost::python::arg("truncation"), boost::python::arg("stepping"),boost::python::arg("newton")), "Set values for the internal configuration of the tracker.  tolerance and truncation are both real doubles.  predictor is a valid value for predictor choice.  stepping and newton are the config structs from pybertini.tracking.config.")

			.def("track_path", &TrackerT::TrackPath, 
				 (boost::python::arg("result"), "start_time", "end_time", "start_point"), 
				 "The main function of the tracker, once its set up.  Feed it, in (result, start_time, end_time, start_point")

			.def("get_system",&TrackerT::GetSystem,return_internal_reference<>(), "Gets an internal reference to the tracked system.")

			.def("predictor",get_predictor_,"Query the current predictor method used by the tracker.")
			.def("predictor",set_predictor_,"Set the predictor method used by the tracker.")

			.def("set_stepsize", &TrackerT::SetStepSize)

			.def("reinitialize_initial_step_size", &TrackerT::ReinitializeInitialStepSize, "Set whether the tracker should re-set the stepsize to the configured-initial stepsize when it starts tracking.  Feed it a bool")
			.def("num_total_steps_taken", &TrackerT::NumTotalStepsTaken,"Ask how many steps have been taken so far, including failures")

			.def("tracking_tolerance", &TrackerT::TrackingTolerance,"A step is labeled as a failure if newton correcting doesn't yield a residual less than this tolerance.  A real number, the smaller the slower tracking, generally speaking")
			.def("tracking_tolerance", &TrackerT::SetTrackingTolerance,"Set the tracking tolerance for the tracker")

			.def("infinite_truncation_tolerance", &TrackerT::SetInfiniteTruncationTolerance,"Set the path truncation tolerance for infinite paths for the tracker")
			.def("infinite_truncation_tolerance", &TrackerT::InfiniteTruncationTolerance,"Get the path truncation tolerance for infinite paths for the tracker")

			.def("infinite_truncation", &TrackerT::SetInfiniteTruncation, "Decide whether the tracker should truncate infinite paths.  See also infinite_truncation_tolerance")
			.def("infinite_truncation", &TrackerT::InfiniteTruncation, "Get the bool for whether the tracker should truncate infinite paths.  See also infinite_truncation_tolerance")

			.def("get_stepping",&TrackerT::template Get<tracking::SteppingConfig>,return_internal_reference<>(),"Get the tracker's internal configuration for things that control stepping behaviour")
			.def("get_newton",&TrackerT::template Get<tracking::NewtonConfig>,return_internal_reference<>(),"Get the tracker's internal configuration for Newton correction")
			.def("set_stepping",&TrackerT::template Set<tracking::SteppingConfig>,"Set the tracker's internal configuration for things that control stepping behaviour")
			.def("set_newton",&TrackerT::template Set<tracking::NewtonConfig>,"Set the tracker's internal configuration for Newton correction")

			.def("current_point", &TrackerT::CurrentPoint)
			.def("current_time", &TrackerT::CurrentTime)
			.def("current_precision", &TrackerT::CurrentPrecision)

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

			.def("refine", return_Refine3_ptr<dbl>)
			.def("refine", return_Refine3_ptr<mpfr>)
			.def("refine", return_Refine4_ptr<dbl>)
			.def("refine", return_Refine4_ptr<mpfr>)
			;
		}

		
		
		template<typename TrackerT>
		template<class PyClass>
		void FixedDoubleTrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("refine", return_Refine3_ptr<dbl>)
			.def("refine", return_Refine4_ptr<dbl>)
			;
		}

		
		template<typename TrackerT>
		template<class PyClass>
		void FixedMultipleTrackerVisitor<TrackerT>::visit(PyClass& cl) const
		{
			cl
			.def("refine", return_Refine3_ptr<mpfr>)
			.def("refine", return_Refine4_ptr<mpfr>)
			;
		}


		template<typename T>
		template<class PyClass>
		void SteppingVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def_readwrite("initial_step_size", &tracking::SteppingConfig::initial_step_size,"The initial stepsize when tracking is started.  See also tracking.AMPTracker.reinitialize_initial_step_size")
			.def_readwrite("max_step_size", &tracking::SteppingConfig::max_step_size,"The maximum allowed stepsize during tracking.  See also min_num_steps")
			.def_readwrite("min_step_size", &tracking::SteppingConfig::min_step_size,"The minimum stepsize the tracker is allowed to take.  See also max_step_size")
			.def_readwrite("step_size_success_factor", &tracking::SteppingConfig::step_size_success_factor,"The scale factor for stepsize, after some consecutive steps.  See also consecutive_successful_steps_before_stepsize_increase")
			.def_readwrite("step_size_fail_factor", &tracking::SteppingConfig::step_size_fail_factor, "The scale factor for stepsize, after a fail happens.  See also step_size_success_factor")
			.def_readwrite("consecutive_successful_steps_before_stepsize_increase", &tracking::SteppingConfig::consecutive_successful_steps_before_stepsize_increase,"This number of successful steps have to taken consecutively, and then the stepsize is permitted to increase")
			.def_readwrite("min_num_steps", &tracking::SteppingConfig::min_num_steps, "The minimum number of steps the tracker can take between now and then.  This is useful if you are tracking closely between times, and want to guarantee some number of steps are taken.  Then again, this could be wasteful, too.")
			.def_readwrite("max_num_steps", &tracking::SteppingConfig::max_num_steps, "The maximum number of steps.  Tracking will die if it tries to take more than this number, sad day.")
			.def_readwrite("frequency_of_CN_estimation", &tracking::SteppingConfig::frequency_of_CN_estimation, "How frequently the condition number should be updated.  Less frequently is faster (estimation requires an additional linear solve), but may cause precision adjustment to lag behind.")
			;
		}

		

		void ExportTrackers()
		{	
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".tracking");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("tracking") = new_submodule;

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Tracking things for PyBertini.  Includes the three fundamental trackers, and utility functions.";

			ExportConfigSettings();
			ExportAMPTracker();
			ExportFixedTrackers();
		}

		void ExportAMPTracker()
		{
			class_<AMPTracker, std::shared_ptr<AMPTracker> >("AMPTracker", "The adaptive multiple precision (AMP) tracker.  Ambient numeric type is multiple-precision (mpfr_complex).  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<AMPTracker>())
			.def(AMPTrackerVisitor<AMPTracker>())
			;
		}

		void ExportFixedTrackers()
		{
			ExportFixedDoubleTracker();
			ExportFixedMultipleTracker();
		}

		void ExportFixedDoubleTracker()
		{
			class_<DoublePrecisionTracker, std::shared_ptr<DoublePrecisionTracker> >("DoublePrecisionTracker", "The double precision tracker.  Tracks using only complex doubles.  Ambient numeric type is double.  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<DoublePrecisionTracker>())
			.def(FixedDoubleTrackerVisitor<DoublePrecisionTracker>())
			;
		}

		void ExportFixedMultipleTracker()
		{
			class_<MultiplePrecisionTracker, std::shared_ptr<MultiplePrecisionTracker> >("MultiplePrecisionTracker", "The fixed multiple precision tracker.  Ambient numeric type is multiple-precision (mpfr_complex).  Precision is the value of pybertini.default_precision() at contruction.  Errors if you try to feed it things not at that precision.  Contruct one by feeding it a system -- cannot be constructed without feeding it a system.  Adjust its settings via configs and the `setup` function.  Then, call method `track_path`.", init<const System&>())
			.def(TrackerVisitor<MultiplePrecisionTracker>())
			.def(FixedMultipleTrackerVisitor<MultiplePrecisionTracker>())
			;

		}
		
		
		
		void ExportConfigSettings()
		{
			using namespace bertini::tracking;

			enum_<Predictor>("Predictor")
				.value("Constant", Predictor::Constant)
				.value("Euler", Predictor::Euler)
				.value("Heun", Predictor::Heun)
				.value("RK4", Predictor::RK4)
				.value("HeunEuler", Predictor::HeunEuler)
				.value("RKNorsett34", Predictor::RKNorsett34)
				.value("RKF45", Predictor::RKF45)
				.value("RKCashKarp45", Predictor::RKCashKarp45)
				.value("RKDormandPrince56", Predictor::RKDormandPrince56)
				.value("RKVerner67", Predictor::RKVerner67)
				;

			enum_<SuccessCode>("SuccessCode")
				.value("Success", SuccessCode::Success)
				.value("HigherPrecisionNecessary", SuccessCode::HigherPrecisionNecessary)
				.value("ReduceStepSize", SuccessCode::ReduceStepSize)
				.value("GoingToInfinity", SuccessCode::GoingToInfinity)
				.value("FailedToConverge", SuccessCode::FailedToConverge)
				.value("MatrixSolveFailure", SuccessCode::MatrixSolveFailure)
				.value("MatrixSolveFailureFirstPartOfPrediction", SuccessCode::MatrixSolveFailureFirstPartOfPrediction)
				.value("MaxNumStepsTaken", SuccessCode::MaxNumStepsTaken)
				.value("MaxPrecisionReached", SuccessCode::MaxPrecisionReached)
				.value("MinStepSizeReached", SuccessCode::MinStepSizeReached)
				.value("Failure", SuccessCode::Failure)
				.value("SingularStartPoint", SuccessCode::SingularStartPoint)
				.value("ExternallyTerminated", SuccessCode::ExternallyTerminated)
				.value("MinTrackTimeReached", SuccessCode::MinTrackTimeReached)
				.value("SecurityMaxNormReached", SuccessCode::SecurityMaxNormReached)
				.value("CycleNumTooHigh", SuccessCode::CycleNumTooHigh)
				;
			
			{ // enter a scope for config types
				scope current_scope;
				std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
				new_submodule_name.append(".config");
				object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
				current_scope.attr("config") = new_submodule;

				scope new_submodule_scope = new_submodule;


				// class_<config::Tolerances<double>>("Tolerances_d",init<>())
				// 	.def(TolerancesVisitor<double>());

				// class_<config::Tolerances<mpfr_float>>("Tolerances_mp",init<>())
				// 	.def(TolerancesVisitor<mpfr_float>());

				class_<SteppingConfig, std::shared_ptr<SteppingConfig> >("SteppingConfig", init<>())
					.def(SteppingVisitor<double>())
					;
				
				
				class_<NewtonConfig, std::shared_ptr<NewtonConfig> >("NewtonConfig", init<>())
					.def_readwrite("max_num_newton_iterations", &NewtonConfig::max_num_newton_iterations)
					.def_readwrite("min_num_newton_iterations", &NewtonConfig::min_num_newton_iterations)
					;
				
				
				class_<FixedPrecisionConfig, std::shared_ptr<FixedPrecisionConfig> >("FixedPrecisionConfig", init<System const&>());


				
				
				class_<AdaptiveMultiplePrecisionConfig, std::shared_ptr<AdaptiveMultiplePrecisionConfig> >("AMPConfig", init<>())
					.def(init<System const&>())
					.def("set_amp_config_from", &AdaptiveMultiplePrecisionConfig::SetAMPConfigFrom)
					.def("set_phi_psi_from_bounds", &AdaptiveMultiplePrecisionConfig::SetPhiPsiFromBounds)
					.def("set_bounds_and_epsilon_from", &AdaptiveMultiplePrecisionConfig::SetBoundsAndEpsilonFrom)
					.def_readwrite("coefficient_bound", &AdaptiveMultiplePrecisionConfig::coefficient_bound)
					.def_readwrite("degree_bound", &AdaptiveMultiplePrecisionConfig::degree_bound)
					.def_readwrite("epsilon", &AdaptiveMultiplePrecisionConfig::epsilon)
					.def_readwrite("phi", &AdaptiveMultiplePrecisionConfig::Phi)
					.def_readwrite("psi", &AdaptiveMultiplePrecisionConfig::Psi)
					.def_readwrite("safety_digits_1", &AdaptiveMultiplePrecisionConfig::safety_digits_1)
					.def_readwrite("safety_digits_2", &AdaptiveMultiplePrecisionConfig::safety_digits_2)
					.def_readwrite("maximum_precision", &AdaptiveMultiplePrecisionConfig::maximum_precision)
					.def_readwrite("consecutive_successful_steps_before_precision_decrease", &AdaptiveMultiplePrecisionConfig::consecutive_successful_steps_before_precision_decrease)
					.def_readwrite("max_num_precision_decreases", &AdaptiveMultiplePrecisionConfig::max_num_precision_decreases)
					.def_readwrite("coefficient_bound", &AdaptiveMultiplePrecisionConfig::coefficient_bound)
					;
				
				def("amp_config_from", &AMPConfigFrom, "make an AMPConfig from a System with generated settings for system-specific things, and default settings otherwise (such as safety digits).");
			}
			
		}



}} // namespaces

