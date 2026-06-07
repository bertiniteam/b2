// python/tracker_config_export.cpp — config structs and enum registration for trackers.
// Separated from tracker_export.cpp so the config registration (cheap) compiles
// independently from the class_<TrackerT> instantiations (expensive).

#include "tracker_export.hpp"

namespace bertini{
	namespace python{

		using namespace bertini::tracking;

		void ExportConfigSettings()
		{
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
				.value("FailedToSelectPrecisionAndStepsize", SuccessCode::FailedToSelectPrecisionAndStepsize)
				;

			{ // enter a scope for config types
				scope current_scope;
				std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
				new_submodule_name.append(".config");
				object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
				current_scope.attr("config") = new_submodule;

				scope new_submodule_scope = new_submodule;

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
