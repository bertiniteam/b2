// python/endgame_config_export.cpp — endgame config struct registration.
// Separated from endgame_{double,mp,amp}_export.cpp so config registration (cheap) compiles
// independently from the class_<EndgameT> instantiations (expensive).

#include "endgame_export.hpp"

namespace bertini{
	namespace python{

		void ExportEndgameSettings()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".config");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("config") = new_submodule;

			scope new_submodule_scope = new_submodule;
			new_submodule_scope.attr("__doc__") = "Endgame configuration structs.";

			class_<endgame::EndgameConfig>("Endgame","Generic endgame settings.  Number of sample points, etc.  Note that some of its configs are rational numbers",init<>())
				.def_readwrite("sample_point_refinement_factor", &endgame::EndgameConfig::sample_point_refinement_factor, "Extra amount of tolerance for refining before computing the final approximation, during endgame.")
				.def_readwrite("num_sample_points", &endgame::EndgameConfig::num_sample_points,"The number of points to use for extrapolant calculation.  In the Power Series Endgame, the is the number of geometrically spaces points on the path.  For Cauchy, this is the number of points on each circle tracked around the target time value.")
				.def_readwrite("min_track_time", &endgame::EndgameConfig::min_track_time,"The minimum distance from the target time to track to.  Decreasing this may help failing runs succeed, or maybe not, because you are, after all, tracking toward a singularity.")
				.def_readwrite("sample_factor", &endgame::EndgameConfig::sample_factor,"The factor by which to space the geometrically spaced 'distance' between sample points, or sample circles for Cauchy.")
				.def_readwrite("max_num_newton_iterations", &endgame::EndgameConfig::max_num_newton_iterations,"the maximum number of newton iterations to be taken during sample point sharpening.  Increasing this can help speed convergence, at the risk of path jumping.")
				.def_readwrite("final_tolerance", &endgame::EndgameConfig::final_tolerance, "The tolerance to which to track the path, using the endgame.  Endgames require two consecutive estimates to be this close to each other under the relative infinity norm.  Default value is 1e-11.")
				;

			class_<endgame::SecurityConfig>("Security","Security settings for endgames.  Control things like truncation because estimated root is near infinity",init<>())
				.def_readwrite("level", &endgame::SecurityConfig::level,"Turns on or off truncation of paths going to infinity during the endgame.  0 is off, 1 is on.")
				.def_readwrite("max_norm", &endgame::SecurityConfig::max_norm,"If on, the norm at which to truncate a path.")
				;

			class_<endgame::PowerSeriesConfig>("PowerSeriesConfig","Settings specific to the power series endgame for computing singular endpoints",init<>())
				.def_readwrite("max_cycle_number", &endgame::PowerSeriesConfig::max_cycle_number,"The maximum cycle number to consider, when calculating the cycle number which best fits the path being tracked.")
				.def_readwrite("cycle_number_amplification", &endgame::PowerSeriesConfig::cycle_number_amplification,"The maximum number allowable iterations during endgames, for points used to approximate the final solution.")
				;

			class_<endgame::CauchyConfig>("CauchyConfig","Settings specific to the Cauchy endgame for computing singular endpoints",init<>())
				.def_readwrite("cycle_cutoff_time", &endgame::CauchyConfig::cycle_cutoff_time)
				.def_readwrite("ratio_cutoff_time", &endgame::CauchyConfig::ratio_cutoff_time)
				.def_readwrite("minimum_for_c_over_k_stabilization", &endgame::CauchyConfig::minimum_for_c_over_k_stabilization)
				.def_readwrite("maximum_cauchy_ratio", &endgame::CauchyConfig::maximum_cauchy_ratio)
				.def_readwrite("num_needed_for_stabilization", &endgame::CauchyConfig::num_needed_for_stabilization,"When running stabilization testing for the cycle number when entering the endgame, this is the number of consecutive points for which the test must pass.")
				.def_readwrite("fail_safe_maximum_cycle_number", &endgame::CauchyConfig::fail_safe_maximum_cycle_number, "max number of loops before giving up." )
				;
		}

}} // namespaces
