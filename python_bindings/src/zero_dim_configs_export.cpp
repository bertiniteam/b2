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
			.def_readwrite("newton_before_endgame", &TolerancesConfig::newton_before_endgame)
			.def_readwrite("newton_during_endgame", &TolerancesConfig::newton_during_endgame)
			.def_readwrite("final_tolerance", &TolerancesConfig::final_tolerance)
			.def_readwrite("path_truncation_threshold", &TolerancesConfig::path_truncation_threshold)
			;

			class_<MidPathConfig>("MidPathConfig", init<>())
			.def_readwrite("same_point_tolerance", &MidPathConfig::same_point_tolerance)
			;

			class_<AutoRetrackConfig>("AutoRetrackConfig", init<>())
			.def_readwrite("midpath_decrease_tolerance_factor", &AutoRetrackConfig::midpath_decrease_tolerance_factor)
			;

			class_<SharpeningConfig>("SharpeningConfig", init<>())
			.def_readwrite("sharpendigits", &SharpeningConfig::sharpendigits)
			.def_readwrite("function_residual_tolerance", &SharpeningConfig::function_residual_tolerance)
			.def_readwrite("ratio_tolerance", &SharpeningConfig::ratio_tolerance)
			;

			class_<RegenerationConfig>("RegenerationConfig", init<>())
			.def_readwrite("remove_infinite_endpoints", &RegenerationConfig::remove_infinite_endpoints)
			.def_readwrite("higher_dimension_check", &RegenerationConfig::higher_dimension_check)
			.def_readwrite("start_level", &RegenerationConfig::start_level)
			.def_readwrite("newton_before_endgame", &RegenerationConfig::newton_before_endgame)
			.def_readwrite("newton_during_endgame", &RegenerationConfig::newton_during_endgame)
			.def_readwrite("final_tolerance", &RegenerationConfig::final_tolerance)
			;

			class_<PostProcessingConfig>("PostProcessingConfig", init<>())
			.def_readwrite("real_threshold", &PostProcessingConfig::real_threshold)
			.def_readwrite("endpoint_finite_threshold", &PostProcessingConfig::endpoint_finite_threshold)
			.def_readwrite("same_point_tolerance", &PostProcessingConfig::same_point_tolerance)
			;

			class_<ZeroDimConfig<dbl_complex>>("ZeroDimConfigDoublePrec", init<>())
			.def_readwrite("start_time", &ZeroDimConfig<dbl_complex>::start_time)
			.def_readwrite("target_time", &ZeroDimConfig<dbl_complex>::target_time)
			.def_readwrite("endgame_boundary", &ZeroDimConfig<dbl_complex>::endgame_boundary)
			;

			class_<ZeroDimConfig<mpfr_complex>>("ZeroDimConfigMultiprec", init<>())
			.def_readwrite("start_time", &ZeroDimConfig<mpfr_complex>::start_time)
			.def_readwrite("target_time", &ZeroDimConfig<mpfr_complex>::target_time)
			.def_readwrite("endgame_boundary", &ZeroDimConfig<mpfr_complex>::endgame_boundary)
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
