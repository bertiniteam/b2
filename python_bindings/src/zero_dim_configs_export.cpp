// python/zero_dim_configs_export.cpp — ZeroDim config structs and metadata registration.
// Separated so the lightweight config/metadata registration compiles independently
// from the heavy class_<ZeroDimT> algorithm instantiations.

#include "zero_dim_export.hpp"
#include <sstream>

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
			.def_readwrite("slice_newton_before_endgame", &RegenerationConfig::slice_newton_before_endgame,
				"Slice-moving tracking tolerance before the endgame (Bertini 1 SliceTolBeforeEG). "
				"Separate from TolerancesConfig.newton_before_endgame, which governs the main tracking.")
			.def_readwrite("slice_newton_during_endgame", &RegenerationConfig::slice_newton_during_endgame,
				"Slice-moving tracking tolerance during the endgame (Bertini 1 SliceTolDuringEG).")
			.def_readwrite("slice_final_tolerance", &RegenerationConfig::slice_final_tolerance,
				"Final tolerance to track the slice move to, using the endgame (Bertini 1 SliceFinalTol).")
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

			// One ZeroDimConfig for every precision model -- the homotopy times are stored precision-free
			// (mpq_rational) and converted to the tracking type at use, so the config is no longer
			// templated on the complex type.  The times are real (the solve tracks the real t-axis);
			// they are exposed as real_mp and round-trip exactly, the same way SteppingConfig exposes
			// its mpq_rational step sizes.
			class_<ZeroDimConfig>("ZeroDimConfig", init<>())
			.add_property("start_time",
				+[](ZeroDimConfig const& c) -> real_mp { return real_mp(c.start_time); },
				+[](ZeroDimConfig& c, real_mp const& v) { c.start_time = mpq_rational(v); },
				"The time value at which the homotopy starts (where the start solutions live).")
			.add_property("target_time",
				+[](ZeroDimConfig const& c) -> real_mp { return real_mp(c.target_time); },
				+[](ZeroDimConfig& c, real_mp const& v) { c.target_time = mpq_rational(v); },
				"The time value the homotopy tracks to (where the solutions of interest live).")
			.add_property("endgame_boundary",
				+[](ZeroDimConfig const& c) -> real_mp { return real_mp(c.endgame_boundary); },
				+[](ZeroDimConfig& c, real_mp const& v) { c.endgame_boundary = mpq_rational(v); },
				"The time value at which tracking stops and the endgame takes over.")
			.def_readwrite("max_num_crossed_path_resolve_attempts", &ZeroDimConfig::max_num_crossed_path_resolve_attempts,
				"How many times to re-track crossed paths (with tightened settings) at the endgame "
				"boundary before giving up. 0 = detect and report only, do not re-track. Default 2.")
			.def_readwrite("num_threads", &ZeroDimConfig::num_threads,
				"Worker threads for a shared-memory (non-MPI) solve. 0 = auto "
				"(all available cores), 1 = serial (no thread pool), N = N threads. The "
				"OMP_NUM_THREADS environment variable overrides this. Threading needs no MPI and "
				"no free-threaded Python: the heavy tracking runs in C++ with the GIL released.")
			;

			// metadata types
			class_<AlgorithmMetaData>("AlgorithmMetaData",init<>())
			.def_readwrite("number_path_failures",&AlgorithmMetaData::number_path_failures)
			.def_readwrite("number_path_successes",&AlgorithmMetaData::number_path_successes)
			.def_readwrite("number_paths_tracked",&AlgorithmMetaData::number_paths_tracked)
			.def_readwrite("start_time",&AlgorithmMetaData::start_time)
			.def_readwrite("elapsed_time",&AlgorithmMetaData::elapsed_time)
			;

			ExposeSolutionMetaData<complex_mp>("SolutionMetaDataMultiPrec");
			ExposeSolutionMetaData<complex_dbl>("SolutionMetaDataDoublePrec");

			ExposeEndgameBoundaryMetaData<complex_mp>("EndgameBoundaryMetaDataMultiPrec");
			ExposeEndgameBoundaryMetaData<complex_dbl>("EndgameBoundaryMetaDataDoublePrec");

			// Report from the path-crossing (midpath) check at the endgame boundary.
			class_<MidpathCheckReport>("MidpathCheckReport", init<>())
			.def_readonly("passed", &MidpathCheckReport::passed,
				"Did the final midpath check pass (no path crossings remained)?  False means one or "
				"more crossings were left unresolved and the affected solutions may be wrong.")
			.def_readonly("num_crossings_detected", &MidpathCheckReport::num_crossings_detected,
				"Number of crossed paths found on the FIRST check, before any re-tracking.")
			.def_readonly("num_resolve_attempts", &MidpathCheckReport::num_resolve_attempts,
				"How many re-track attempts were actually performed.")
			.add_property("crossed_path_indices",
				+[](MidpathCheckReport const& r){
					boost::python::list out;
					for (auto i : r.crossed_path_indices) out.append(i);
					return out;
				},
				"Indices of the paths flagged as crossed on the first check.")
			;

			// SolveReport: the end-of-solve diagnostic summary (see the solver's report() method).
			class_<SolveReport>("SolveReport", init<>())
			.def_readonly("num_paths_tracked", &SolveReport::num_paths_tracked,
				"Total number of paths tracked (the start-system / Bezout count).")
			.def_readonly("num_finite_solutions", &SolveReport::num_finite_solutions,
				"Number of DISTINCT finite solutions (multiple roots counted once).")
			.def_readonly("num_finite_endpoints", &SolveReport::num_finite_endpoints,
				"Raw number of finite, successful endpoints (before collapsing multiplicities).")
			.def_readonly("num_diverged", &SolveReport::num_diverged,
				"Number of paths that diverged to infinity -- a result, not a failure.")
			.def_readonly("num_failed", &SolveReport::num_failed,
				"Number of paths the tracker could not resolve -- each one a possibly-missing solution.")
			.def_readonly("num_singular", &SolveReport::num_singular,
				"Number of finite solutions flagged singular (multiple or ill-conditioned).")
			.def_readonly("num_real", &SolveReport::num_real,
				"Number of finite solutions flagged real.")
			.def_readonly("max_condition_number", &SolveReport::max_condition_number,
				"Largest condition number among the finite solutions.")
			.def_readonly("max_precision_used", &SolveReport::max_precision_used,
				"Highest working precision (digits) any path needed.")
			.def_readonly("midpath", &SolveReport::midpath,
				"The MidpathCheckReport from the path-crossing check.")
			.def_readonly("all_paths_resolved", &SolveReport::all_paths_resolved,
				"True iff no path failed and no crossing was left unresolved -- the solve is trustworthy.")
			.add_property("failures_by_reason",
				+[](SolveReport const& r){
					boost::python::dict out;
					for (auto const& kv : r.failures_by_reason) out[kv.first] = kv.second;
					return out;
				},
				"Dict {SuccessCode: count} of how the failed paths ended.")
			.def("__str__",  +[](SolveReport const& r){ std::ostringstream s; s << r; return s.str(); })
			.def("__repr__", +[](SolveReport const& r){ std::ostringstream s; s << r; return s.str(); })
			;
		}

}} // namespaces
