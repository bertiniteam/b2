// python/endgame_config_export.cpp — endgame config struct registration.
// Separated from endgame_{double,mp,amp}_export.cpp so config registration (cheap) compiles
// independently from the class_<EndgameT> instantiations (expensive).

#include "endgame_export.hpp"

namespace bertini{
    namespace python{

        void ExportEndgameSettings()
        {
            // config classes live directly in the endgame module, with Config-suffixed names
            // (the .config submodule was flattened away)
            class_<endgame::EndgameConfig>("EndgameConfig","Generic endgame settings.  Number of sample points, etc.  Note that some of its configs are rational numbers",init<>())
                .def_readwrite("sample_point_refinement_factor", &endgame::EndgameConfig::sample_point_refinement_factor, "Extra amount of tolerance for refining before computing the final approximation, during endgame.")
                .def_readwrite("num_sample_points", &endgame::EndgameConfig::num_sample_points,"The number of points to use for extrapolant calculation.  In the Power Series Endgame, the is the number of geometrically spaces points on the path.  For Cauchy, this is the number of points on each circle tracked around the target time value.")
                .def_readwrite("min_track_time", &endgame::EndgameConfig::min_track_time,"The minimum distance from the target time to track to.  Decreasing this may help failing runs succeed, or maybe not, because you are, after all, tracking toward a singularity.")
                .def_readwrite("sample_factor", &endgame::EndgameConfig::sample_factor,"The factor by which to space the geometrically spaced 'distance' between sample points, or sample circles for Cauchy.")
                .def_readwrite("max_num_refinements", &endgame::EndgameConfig::max_num_refinements,"the maximum number of Newton refinements to be taken during sample point sharpening.  Increasing this can help speed convergence, at the risk of path jumping.")
                .def_readwrite("final_tolerance", &endgame::EndgameConfig::final_tolerance, "The tolerance to which to track the path, using the endgame.  Endgames require two consecutive estimates to be this close to each other under the relative infinity norm.  Default value is 1e-11.")
                .def_readwrite("refine_when_increasing_precision", &endgame::EndgameConfig::refine_when_increasing_precision,
                    "Whether to re-refine the samples the endgame is keeping when it moves to a higher "
                    "working precision (default False).  Off, the retained samples carry the accuracy they "
                    "were computed at, which is the accuracy of the lower precision; on, they are refined "
                    "again at the new one, which costs Newton steps per retained sample and buys a sharper "
                    "sample window.  Part of the run's identity, so changing it is a different ask.")
                .def_readwrite("minimum_for_c_over_k_stabilization", &endgame::EndgameConfig::minimum_for_c_over_k_stabilization,
                    "How closely successive c/k estimates must agree for the path to count as being in the "
                    "endgame operating zone, as a ratio of the smaller to the larger.  Their settling is what "
                    "says the cycle-number estimate has settled, and so that the Puiseux asymptotics the "
                    "endgame is built on dominate.")
                .def_readwrite("num_needed_for_stabilization", &endgame::EndgameConfig::num_needed_for_stabilization,
                    "How many consecutive c/k estimates must agree for the path to count as being in the "
                    "endgame operating zone.")
                ;

            class_<endgame::SecurityConfig>("SecurityConfig","Security settings for endgames.  Control things like truncation because estimated root is near infinity",init<>())
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
                .def_readwrite("maximum_cauchy_ratio", &endgame::CauchyConfig::maximum_cauchy_ratio)
                .def_readwrite("fail_safe_maximum_cycle_number", &endgame::CauchyConfig::fail_safe_maximum_cycle_number, "max number of loops before giving up." )
                .def_readwrite("num_consecutive_same_cycle_number", &endgame::CauchyConfig::num_consecutive_same_cycle_number, "Number of consecutive Cauchy approximations that must report the same cycle number before a converged approximation is trusted." )
                ;
        }

}} // namespaces
