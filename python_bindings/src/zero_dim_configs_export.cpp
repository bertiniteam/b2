// python/zero_dim_configs_export.cpp — ZeroDim config structs and metadata registration.
// Separated so the lightweight config/metadata registration compiles independently
// from the heavy class_<ZeroDimT> algorithm instantiations.

#include "zero_dim_export.hpp"

namespace bertini{
	namespace python{

		void ExportZDConfigs()
		{
			using namespace bertini::algorithm;

			class_<TolerancesConfig>("TolerancesConfig", init<>())
			.def_readwrite("newton_before_endgame", &TolerancesConfig::newton_before_endgame,
				"Tracking (Newton) tolerance used while tracking before the endgame begins. "
				"Tighten this if paths drift together or are missed.")
			.def_readwrite("newton_during_endgame", &TolerancesConfig::newton_during_endgame,
				"Tracking (Newton) tolerance used during the endgame.")
			.def_readwrite("final_tolerance", &TolerancesConfig::final_tolerance,
				"The tolerance to which a solution is computed by the endgame. The same-point test "
				"derives its tolerance from this (see PostProcessingConfig.same_point_tolerance_multiplier).")
			.def_readwrite("path_truncation_threshold", &TolerancesConfig::path_truncation_threshold,
				"If a path point's norm exceeds this, the tracker declares the path divergent and stops it.")
			;

			class_<MidPathConfig>("MidPathConfig", init<>())
			.def_readwrite("same_point_tolerance", &MidPathConfig::same_point_tolerance,
				"Tolerance used by the midpath check to detect two paths that have crossed "
				"(become the same point) partway through tracking, which signals a need to retrack.")
			;

			class_<AutoRetrackConfig>("AutoRetrackConfig", init<>())
			.def_readwrite("midpath_decrease_tolerance_factor", &AutoRetrackConfig::midpath_decrease_tolerance_factor,
				"Factor by which tracking tolerances are tightened when retracking after the midpath "
				"check detects a path crossing.")
			;

			class_<SharpeningConfig>("SharpeningConfig", init<>())
			.def_readwrite("sharpendigits", &SharpeningConfig::sharpendigits,
				"How many digits should be correct after sharpening a solution.")
			.def_readwrite("function_residual_tolerance", &SharpeningConfig::function_residual_tolerance,
				"A function value is considered zero if its magnitude is smaller than this.")
			.def_readwrite("ratio_tolerance", &SharpeningConfig::ratio_tolerance,
				"A value is considered zero if the ratio of two different approximations is smaller than this.")
			;

			class_<RegenerationConfig>("RegenerationConfig", init<>())
			.def_readwrite("remove_infinite_endpoints", &RegenerationConfig::remove_infinite_endpoints,
				"Whether endpoints found to be at infinity during regeneration start-point buildup "
				"are discarded. Set True if you are not interested in solutions at infinity.")
			.def_readwrite("higher_dimension_check", &RegenerationConfig::higher_dimension_check,
				"Whether to test for, and remove, points lying on higher-dimensional components during regeneration.")
			.def_readwrite("start_level", &RegenerationConfig::start_level,
				"The regeneration level at which to begin.")
			.def_readwrite("newton_before_endgame", &RegenerationConfig::newton_before_endgame,
				"Regeneration slice tracking tolerance before the endgame.")
			.def_readwrite("newton_during_endgame", &RegenerationConfig::newton_during_endgame,
				"Regeneration slice tracking tolerance during the endgame.")
			.def_readwrite("final_tolerance", &RegenerationConfig::final_tolerance,
				"Regeneration slice final tolerance, tracked to using the endgame.")
			;

			class_<PostProcessingConfig>("PostProcessingConfig", init<>())
			.def_readwrite("real_threshold", &PostProcessingConfig::real_threshold,
				"Bertini 1's ImagThreshold. A (dehomogenized) endpoint is classified real if the "
				"infinity norm of its coordinates' imaginary parts is below this. Default 1e-8.")
			.def_readwrite("endpoint_finite_threshold", &PostProcessingConfig::endpoint_finite_threshold,
				"Bertini 1's EndpointFiniteThreshold. An endpoint is classified at infinity if the "
				"infinity norm of its dehomogenized coordinates exceeds this value. Default 1e5.")
			.def_readwrite("same_point_tolerance_multiplier", &PostProcessingConfig::same_point_tolerance_multiplier,
				"Bertini 1's EndpointSameThreshold. A multiplier (>= 1) on final_tolerance: two "
				"endpoints are the same point (raising multiplicity) when the infinity norm of the "
				"difference of their dehomogenized coordinates is below "
				"final_tolerance * same_point_tolerance_multiplier. Default 10.")
			.def_readwrite("condition_number_threshold", &PostProcessingConfig::condition_number_threshold,
				"Bertini 1's CondNumThreshold. An endpoint is classified singular if it is the endpoint "
				"of multiple paths (multiplicity > 1), or if its spectral-norm condition-number estimate "
				"exceeds this value. Default 1e8.")
			;

			class_<ZeroDimConfig<dbl_complex>>("ZeroDimConfigDoublePrec", init<>())
			.def_readwrite("start_time", &ZeroDimConfig<dbl_complex>::start_time,
				"The time value at which the homotopy starts (where the start solutions live).")
			.def_readwrite("target_time", &ZeroDimConfig<dbl_complex>::target_time,
				"The time value the homotopy tracks to (where the solutions of interest live).")
			.def_readwrite("endgame_boundary", &ZeroDimConfig<dbl_complex>::endgame_boundary,
				"The time value at which tracking stops and the endgame takes over.")
			;

			class_<ZeroDimConfig<mpfr_complex>>("ZeroDimConfigMultiprec", init<>())
			.def_readwrite("start_time", &ZeroDimConfig<mpfr_complex>::start_time,
				"The time value at which the homotopy starts (where the start solutions live).")
			.def_readwrite("target_time", &ZeroDimConfig<mpfr_complex>::target_time,
				"The time value the homotopy tracks to (where the solutions of interest live).")
			.def_readwrite("endgame_boundary", &ZeroDimConfig<mpfr_complex>::endgame_boundary,
				"The time value at which tracking stops and the endgame takes over.")
			;

			// metadata types
			class_<AlgorithmMetaData>("AlgorithmMetaData",init<>())
			.def_readwrite("number_path_failures",&AlgorithmMetaData::number_path_failures)
			.def_readwrite("number_path_successes",&AlgorithmMetaData::number_path_successes)
			.def_readwrite("number_paths_tracked",&AlgorithmMetaData::number_paths_tracked)
			.def_readwrite("start_time",&AlgorithmMetaData::start_time)
			.def_readwrite("elapsed_time",&AlgorithmMetaData::elapsed_time)
			;

			ExposeSolutionMetaData<mpfr_complex>("SolutionMetaDataMultiPrec");
			ExposeSolutionMetaData<dbl_complex>("SolutionMetaDataDoublePrec");

			ExposeEndgameBoundaryMetaData<mpfr_complex>("EndgameBoundaryMetaDataMultiPrec");
			ExposeEndgameBoundaryMetaData<dbl_complex>("EndgameBoundaryMetaDataDoublePrec");
		}

}} // namespaces
