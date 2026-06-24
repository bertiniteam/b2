//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/zero_dim_solve.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/zero_dim_solve.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/nag_algorithms/zero_dim_solve.hpp

\brief Provides the algorithm for computing all zero-dimensional solutions for an algberaic system.
*/

#pragma once

#include "bertini2/num_traits.hpp"

#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"
#include "bertini2/nag_algorithms/midpath_check.hpp"
#include "bertini2/io/generators.hpp"

#include "bertini2/detail/configured.hpp"
#include "bertini2/detail/observable.hpp"

#include "bertini2/nag_algorithms/common/algorithm_base.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"
#include "bertini2/nag_algorithms/events.hpp"
#include "bertini2/nag_algorithms/common/policies.hpp"
#include "bertini2/parallel.hpp"
#include <chrono>
#include <mutex>
#include <iostream>
#include <map>


namespace bertini {

	namespace algorithm {


/**
forward declare of ZeroDim algorithm
*/
template<	typename TrackerType, typename EndgameType,
			typename SystemType, typename StartSystemType,
			template<typename,typename> class SystemManagementP = policy::CloneGiven >
struct ZeroDim;



/**
specify the traits for the algorithm.  this is why we need the forward declare
*/
template<typename TrackerType, typename EndgameType,
			typename SystemType, typename StartSystemType,
			template<typename,typename> class SystemManagementP>
struct AlgoTraits <ZeroDim<TrackerType, EndgameType, SystemType, StartSystemType, SystemManagementP>>
{
	using BaseRealT = typename tracking::TrackerTraits<TrackerType>::BaseRealT;
	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;

	using NeededConfigs = detail::TypeList<
								TolerancesConfig,
								PostProcessingConfig,
								ZeroDimConfig,
								AutoRetrackConfig
								>;
};


struct AnyZeroDim : public virtual AnyAlgorithm
{
	virtual void WriteMainData(std::ostream& out) const = 0;
	virtual void WriteRawData(std::ostream& out)  const = 0;
	virtual void ApplyParsedConfigs(std::string const& config_str) = 0;
	virtual ~AnyZeroDim() = default;
};




using SolnIndT = typename SolnCont<dbl_complex>::size_type;


/// metadata structs

struct AlgorithmMetaData
{	
	
	SolnIndT number_path_failures = 0;
	SolnIndT number_path_successes = 0;
	SolnIndT number_paths_tracked = 0;

	std::chrono::system_clock::time_point start_time;
	std::chrono::microseconds elapsed_time;
};


template<typename ComplexT>
struct SolutionMetaData
{
	using SolnIndT = typename SolnCont<ComplexT>::size_type;

	// only vaguely metadata.  artifacts of randomness or ordering
	SolnIndT path_index;     		// path number of the solution
	SolnIndT solution_index;      	// solution number

	///// things computed across all of the solve
	bool precision_changed = false;
	ComplexT time_of_first_prec_increase;    // time value of the first increase in precision
	decltype(DefaultPrecision()) max_precision_used = 0;

	///// things computed in pre-endgame only
	SuccessCode pre_endgame_success = SuccessCode::NeverStarted;     // success code


	///// things computed in endgame only
	NumErrorT condition_number; 				// the latest estimate on the condition number
	NumErrorT newton_residual; 				// the latest newton residual
	ComplexT final_time_used;   			// the final value of time tracked to
	NumErrorT accuracy_estimate; 			// accuracy estimate between extrapolations
	NumErrorT accuracy_estimate_user_coords;	// accuracy estimate between extrapolations, in natural coordinates
	unsigned cycle_num;    						// cycle number used in extrapolations
	SuccessCode endgame_success = SuccessCode::NeverStarted;      // success code


	///// things added by post-processing
	NumErrorT function_residual; 	// the latest function residual

	int multiplicity = 1; 		// multiplicity
	bool is_real = false;       		// real flag: whether the (dehomogenized) endpoint is real
	bool is_finite = false;     		// finite flag: whether the endpoint is finite (not at infinity)
	bool is_singular = false;       		// singular flag: whether the endpoint is singular (multiple, or ill-conditioned)

	bool operator==(const SolutionMetaData<ComplexT> & other){ 
		bool result = 
			this->path_index == other.path_index
			 && this->solution_index == other.solution_index
			 && this->precision_changed == other.precision_changed
			 && this->time_of_first_prec_increase == other.time_of_first_prec_increase
			 && this->max_precision_used == other.max_precision_used
			 && this->pre_endgame_success == other.pre_endgame_success
			 && this->condition_number == other.condition_number
			 && this->newton_residual == other.newton_residual
			 && this->final_time_used == other.final_time_used
			 && this->accuracy_estimate == other.accuracy_estimate
			 && this->accuracy_estimate_user_coords == other.accuracy_estimate_user_coords
			 && this->cycle_num == other.cycle_num
			 && this->endgame_success == other.endgame_success
			 && this->function_residual == other.function_residual
			 && this->multiplicity == other.multiplicity
			 && this->is_real == other.is_real
			 && this->is_finite == other.is_finite
			 && this->is_singular == other.is_singular
		;

		return result; }
};

// this is for interoperability with vectors of these in the Python bindings, for better or for worse.
template<typename NumT>
std::ostream& operator<<(std::ostream & out, const SolutionMetaData<NumT> & meta){
	out << "path_index = " << meta.path_index << std::endl;
	out << "solution_index = " << meta.solution_index << std::endl;

	out << "precision_changed = " << meta.precision_changed << std::endl;
	out << "time_of_first_prec_increase = " << meta.time_of_first_prec_increase << std::endl;
	out << "max_precision_used = " << meta.max_precision_used << std::endl;

	out << "pre_endgame_success = " << meta.pre_endgame_success << std::endl;

	out << "condition_number = " << meta.condition_number << std::endl;
	out << "newton_residual = " << meta.newton_residual << std::endl;
	out << "final_time_used = " << meta.final_time_used << std::endl;
	out << "accuracy_estimate = " << meta.accuracy_estimate << std::endl;
	out << "accuracy_estimate_user_coords = " << meta.accuracy_estimate_user_coords << std::endl;
	out << "cycle_num = " << meta.cycle_num << std::endl;
	out << "endgame_success = " << meta.endgame_success << std::endl;

	out << "function_residual = " << meta.function_residual << std::endl;

	out << "multiplicity = " << meta.multiplicity << std::endl;
	out << "is_real = " << meta.is_real << std::endl;
	out << "is_finite = " << meta.is_finite << std::endl;
	out << "is_singular = " << meta.is_singular << std::endl;

	return out;
}


template<typename ComplexT>
struct EGBoundaryMetaData
{	
	using RealT = typename NumTraits<ComplexT>::Real;

	Vec<ComplexT> path_point;
	SuccessCode success_code = SuccessCode::NeverStarted;
	RealT last_used_stepsize;
	// The precision the tracker was actually using when it reached the endgame boundary.
	// Carried explicitly (rather than inferred from Precision(path_point)) because the
	// path_point is always stored as the tracker's BaseComplexT (multiprecision for AMP);
	// a path that tracked in double gets widened on output, so its mantissa precision no
	// longer reflects the precision that was in use.  The endgame should resume at this
	// precision.  See zero_dim_solve TrackSinglePathDuringEG.
	unsigned precision = DoublePrecision();

	EGBoundaryMetaData() = default;
	EGBoundaryMetaData(EGBoundaryMetaData const&) = default;
	EGBoundaryMetaData(Vec<ComplexT> const& pt, SuccessCode const& code, RealT const& ss, unsigned prec) :
		path_point(pt), success_code(code), last_used_stepsize(ss), precision(prec)
	{}

	bool operator==(const EGBoundaryMetaData<ComplexT> & other){
		bool result =
			this->path_point == other.path_point
			&& this->success_code == other.success_code
			&& this->last_used_stepsize == other.last_used_stepsize
			&& this->precision == other.precision
		;

		return result;
	}
};

// this is for interoperability with vectors of these in the Python bindings, for better or for worse.
template<typename NumT>
std::ostream& operator<<(std::ostream & out, const EGBoundaryMetaData<NumT> & meta){
	out << "path_point = " << meta.path_point << std::endl;
	out << "success_code = " << meta.success_code << std::endl;
	out << "last_used_stepsize = " << meta.last_used_stepsize << std::endl;
	out << "precision = " << meta.precision << std::endl;
	return out;
}

/**
\brief Summary of what the midpath (path-crossing) check found at the endgame boundary, and what the
algorithm did about it.

Two distinct paths landing on the same point at the endgame boundary is a probability-0 event: it
signals a path crossing (under-resolved tracking), not a benign coincidence.  The zero-dim algorithm
detects these at the boundary using a relaxed same-point tolerance and re-tracks the offending paths
with tightened settings (see ZeroDim::EGBoundaryAction).  This struct records the outcome so a caller
can tell whether a solve hit crossings and whether they were resolved.

\see ZeroDim::EndgameBoundaryMetadata
*/
struct MidpathCheckReport
{
	bool passed = true;                                    ///< Did the final midpath check pass (no crossings remained)?
	unsigned num_crossings_detected = 0;                   ///< Number of crossed paths found on the *first* check, before any re-tracking.
	unsigned num_resolve_attempts = 0;                     ///< How many re-track attempts were actually performed.
	std::vector<unsigned long long> crossed_path_indices;  ///< Indices of the paths flagged as crossed on the first check.
};

inline
std::ostream& operator<<(std::ostream & out, const MidpathCheckReport & r){
	out << "passed = " << std::boolalpha << r.passed << std::endl;
	out << "num_crossings_detected = " << r.num_crossings_detected << std::endl;
	out << "num_resolve_attempts = " << r.num_resolve_attempts << std::endl;
	out << "crossed_path_indices = [";
	for (size_t i = 0; i < r.crossed_path_indices.size(); ++i)
		out << (i ? ", " : "") << r.crossed_path_indices[i];
	out << "]" << std::endl;
	return out;
}

/**
\brief A concise, human-readable summary of a zero-dimensional solve: how every path ended up, and
-- crucially -- whether any path was lost.

A raw solution count can lie: if a near-singular path fails to track, a genuine root is silently
missing and the count comes back short.  This report buckets every endpoint into a finite solution,
a clean divergence, or a genuine FAILURE (the tracker gave up), records the failures by their named
SuccessCode, and exposes all_paths_resolved as the one-line health check.

\see ZeroDim::Report, SummarizeSolve
*/
struct SolveReport
{
	unsigned long long num_paths_tracked = 0;    ///< total paths tracked (the Bezout / start count)
	unsigned long long num_finite_solutions = 0; ///< DISTINCT finite solutions (round of sum 1/multiplicity)
	unsigned long long num_finite_endpoints = 0; ///< raw count of finite, successful endpoints
	unsigned long long num_diverged = 0;         ///< paths ending at infinity (Success & !finite, or GoingToInfinity)
	unsigned long long num_failed = 0;           ///< paths the tracker could not resolve (no solution, no clean divergence)
	unsigned long long num_singular = 0;         ///< finite solutions flagged singular (multiple / ill-conditioned)
	unsigned long long num_real = 0;             ///< finite solutions flagged real
	std::map<SuccessCode, unsigned long long> failures_by_reason; ///< histogram of the failed paths' SuccessCodes
	double max_condition_number = 0;             ///< largest condition number among finite solutions
	unsigned max_precision_used = 0;             ///< highest working precision any path needed (digits)
	MidpathCheckReport midpath;                  ///< the path-crossing check outcome
	bool all_paths_resolved = false;             ///< num_failed==0 AND no unresolved crossings: the solve is trustworthy
};

/**
\brief Build a SolveReport from the per-endpoint final metadata and the midpath-crossing report.
Pure aggregation -- one pass, no solve state -- so it is unit-testable on synthetic metadata.
*/
template<typename ComplexT>
SolveReport SummarizeSolve(std::vector<SolutionMetaData<ComplexT>> const& metadata,
                           MidpathCheckReport const& midpath)
{
	SolveReport r;
	r.num_paths_tracked = metadata.size();
	r.midpath = midpath;
	double finite_distinct = 0;
	for (auto const& m : metadata)
	{
		if (m.max_precision_used > r.max_precision_used)
			r.max_precision_used = m.max_precision_used;

		bool diverged = (m.endgame_success == SuccessCode::GoingToInfinity);
		if (m.endgame_success == SuccessCode::Success)
		{
			if (m.is_finite)
			{
				++r.num_finite_endpoints;
				finite_distinct += 1.0 / m.multiplicity;
				if (m.is_singular) ++r.num_singular;
				if (m.is_real)     ++r.num_real;
				if (static_cast<double>(m.condition_number) > r.max_condition_number)
					r.max_condition_number = static_cast<double>(m.condition_number);
			}
			else
				diverged = true;
		}

		if (diverged)
			++r.num_diverged;
		else if (m.endgame_success != SuccessCode::Success)   // neither a solution nor a clean divergence
		{
			++r.num_failed;
			++r.failures_by_reason[m.endgame_success];
		}
	}
	r.num_finite_solutions = static_cast<unsigned long long>(finite_distinct + 0.5);
	r.all_paths_resolved = (r.num_failed == 0) && midpath.passed;
	return r;
}

inline
std::ostream& operator<<(std::ostream & out, const SolveReport & r)
{
	out << "zero-dim solve -- " << r.num_paths_tracked << " paths tracked\n";
	out << "  finite solutions    " << r.num_finite_solutions << "   (distinct)\n";
	out << "  diverged            " << r.num_diverged << "\n";
	out << "  FAILED              " << r.num_failed;
	if (r.num_failed)
	{
		out << "    ";
		bool first = true;
		for (auto const& kv : r.failures_by_reason)
		{
			out << (first ? "" : ", ") << kv.second << " x " << kv.first;
			first = false;
		}
	}
	out << "\n  ----\n";
	out << "  singular solutions  " << r.num_singular << "\n";
	out << "  real solutions      " << r.num_real << "\n";
	out << "  path crossings      " << r.midpath.num_crossings_detected
	    << (r.midpath.passed ? " (resolved)" : " (UNRESOLVED)") << "\n";
	out << "  max condition num   " << r.max_condition_number << "\n";
	out << "  max precision used  " << r.max_precision_used << " digits\n";
	out << "  all paths resolved? " << (r.all_paths_resolved ? "yes" : "NO") << "\n";
	return out;
}

/**
\brief the basic zero dim algorithm, which solves a system.
*/
		template<	typename TrackerType, typename EndgameType,
					typename SystemType, typename StartSystemType,
					template<typename,typename> class SystemManagementP>
		struct ZeroDim :
							public virtual AnyZeroDim,
							public Observable,
							public SystemManagementP<SystemType, StartSystemType>,
							public detail::Configured<
								typename AlgoTraits< ZeroDim<TrackerType, EndgameType, SystemType, StartSystemType, SystemManagementP>>::NeededConfigs>
		{
			// these usings are for getters in python
			using TrackerT          = TrackerType;
			using EndgameT          = EndgameType;
			using SystemT           = SystemType;
			using StartSystemT       = StartSystemType;

			// This algorithm emits its lifecycle events on the AnyZeroDim base, so a
			// single observer type can watch any templated ZeroDim.  Accept observers
			// declared for AnyZeroDim (in addition to the exact concrete type).
			bool ObservableIsA(std::type_index t) const override
			{
				return Observable::ObservableIsA(t) || t == std::type_index(typeid(AnyZeroDim));
			}




/// a bunch of using statements to reduce typing.
			using BaseComplexT 	= typename tracking::TrackerTraits<TrackerType>::BaseComplexT;
			using BaseRealT    	= typename tracking::TrackerTraits<TrackerType>::BaseRealT;

			using PrecisionConfig 	= typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

			using SolnIndT 			= typename SolnCont<BaseComplexT>::size_type;

			using SystemManagementPolicy = SystemManagementP<SystemType, StartSystemType>;

			using StoredSystemT = typename SystemManagementPolicy::StoredSystemT;
			using StoredStartSystemT = typename SystemManagementPolicy::StoredStartSystemT;


			using Config = detail::Configured<
								typename AlgoTraits<ZeroDim<TrackerType, EndgameType, SystemType, StartSystemType, SystemManagementP>>::NeededConfigs>;
			using Config::Get;


			using Tolerances = TolerancesConfig;
			using PostProcessing = PostProcessingConfig;
			using ZeroDimConf = ZeroDimConfig;
			using AutoRetrack = AutoRetrackConfig;


			using EGBoundaryMetaDataT = EGBoundaryMetaData<BaseComplexT>;
			using SolutionMetaDataT = SolutionMetaData<BaseComplexT>;


// a few more using statements

			using MidpathType = MidpathChecker<BaseRealT, BaseComplexT, EGBoundaryMetaData<BaseComplexT>>;

			using SystemManagementPolicy::TargetSystem;
			using SystemManagementPolicy::StartSystem;
			using SystemManagementPolicy::Homotopy;



/// constructors

			/**
			Construct a ZeroDim algorithm object.

			You must at least pass in the system used to track, though the particular arguments required depend on the policies used in your instantiation of ZeroDim.

			\see RefToGiven, CloneGiven
			*/
			template<typename ... SysTs>
			ZeroDim(SysTs const& ...sys) : SystemManagementPolicy(sys...), tracker_(TargetSystem()), endgame_(tracker_)
			{
				ConsistencyCheck();
				DefaultSetup();
			}


			/**
			\brief Main Run() function provided for calling from the blackbox mode
			*/
			void Run() override
			{
				// TODO(Phase 2 Python): AnyZeroDim lacks GetSolutions() -- solutions are only
				// accessible via the concrete type's member arrays or WriteMainData().
				// Resolve when designing Python bindings; the right approach depends on whether
				// Python holds AnyZeroDim or the concrete ZeroDim<...> type directly.
#ifdef BERTINI2_HAVE_MPI
				if (parallel::Size() > 1)
				{
					RunParallel(parallel::WorldComm());
					return;
				}
#endif
				Solve();
			}

#ifdef BERTINI2_HAVE_MPI
			/**
			\brief Manager-worker parallel solve using plain C MPI.

			All ranks must have already constructed an identical ZeroDim from the same input.
			Rank 0 acts as manager; ranks 1..P-1 are workers.

			Phase 1: workers track paths to the endgame boundary.
			         Manager gathers results, runs EGBoundaryAction (midpath check).
			Phase 2: manager distributes successful paths (with boundary data) to workers.
			         Workers run the endgame, send final solutions back.
			Manager runs PostEGAction (multiplicities, same-point detection).
			*/
			void RunParallel(MPI_Comm comm)
			{
				using Result = parallel::FullPathResult<BaseComplexT>;
				using Task   = parallel::StartPointTask<BaseComplexT>;

				mpfr_free_cache(); // reproducibility: clear this rank's mpfr constant cache (see Solve())

				solutions_user_coords_fresh_ = false;

				// Each rank built its homotopy independently in the constructor -- with its OWN
				// random patch, start-system coefficients, and gamma -- so a given path index would
				// mean a different path on each worker, and the manager would gather a scrambled,
				// mostly-wrong solution set.  Make rank 0 the single authoritative source: broadcast
				// its actual target system, start system, and homotopy to every rank.  These are all
				// exact (rational coefficients + integer/rational patch), so serialization carries
				// them bit-for-bit -- there is no per-rank re-derivation that could drift, and the
				// distributed homotopy is identical to rank 0's serial one.  We also broadcast rank
				// 0's RNG seed so any randomness drawn later (e.g. in the endgame) is consistent; the
				// per-path tracking RNG is already deterministic in the path index (ReseedThisThread).
				{
					unsigned long seed = parallel::IsManager() ? GetGlobalSeed() : 0ul;
					MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, 0, comm);
					SetGlobalSeed(seed);

					// When this policy owns its systems (CloneGiven), install rank 0's authoritative
					// systems on every rank.  For RefToGiven the user manages the systems and is
					// responsible for their consistency across ranks, so we leave them untouched
					// (matching the old no-op SystemSetup for that policy).
					if constexpr (SystemManagementPolicy::OwnsSystems)
					{
						parallel::mpi_broadcast_serialized(comm, TargetSystem(), 0);
						parallel::mpi_broadcast_serialized(comm, StartSystem(), 0);
						parallel::mpi_broadcast_serialized(comm, Homotopy(),    0);
					}

					num_start_points_ = StartSystem().NumStartPoints();
					GetTracker().SetSystem(Homotopy());
				}

				PreSolveChecks();
				PreSolveSetup();

				const int n_threads = parallel::WorkerThreadCount();

				// One round = dispatch a set of path indices, each executed as a WHOLE path
				// (pre-endgame + endgame) by a worker.  Round 0 is every path; later rounds re-run
				// only the crossed paths with escalated settings -- the re-track is path-independent,
				// so it rides the same demand-driven worker pool (no idle-manager bottleneck).
				if (parallel::IsManager())
				{
					// Rank 0 computes every start point once, authoritatively, and ships each to the
					// worker that tracks it.  Cached so re-track rounds reuse the identical points.
					std::vector<Vec<BaseComplexT>> start_points(num_start_points_);
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
						start_points[ii] = ComputeStartPoint(static_cast<SolnIndT>(ii));

					auto run_round = [&](std::queue<Task>& q)
					{
						parallel::RunManagerLoop<Task, Result>(comm, q,
							[this](Result const& r){ StoreFullPathResult(r); });
					};

					std::queue<Task> queue;
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
						queue.push(Task{ static_cast<SolnIndT>(ii), start_points[ii] });
					run_round(queue);

					// Midpath/crossing check on the collected boundary points, then bounded parallel
					// resolve rounds.  A continue/stop flag is broadcast to the workers each round so
					// they know whether to run again (and escalate their settings first).
					auto passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
					midpath_report_ = MidpathCheckReport{};
					midpath_report_.passed = passed;
					for (auto const& v : midpath_.GetCrossedPaths())
						midpath_report_.crossed_path_indices.push_back(v.index());
					midpath_report_.num_crossings_detected =
						static_cast<unsigned>(midpath_report_.crossed_path_indices.size());

					const auto max_attempts = this->template Get<ZeroDimConf>().max_num_crossed_path_resolve_attempts;
					unsigned num_resolve_attempts = 0;
					while (true)
					{
						int keep_going = (!passed && num_resolve_attempts < max_attempts) ? 1 : 0;
						MPI_Bcast(&keep_going, 1, MPI_INT, 0, comm);
						if (!keep_going)
							break;

						std::queue<Task> redo;
						for (auto const& v : midpath_.GetCrossedPaths())
							if (v.rerun())
							{
								auto idx = static_cast<SolnIndT>(v.index());
								redo.push(Task{ idx, start_points[idx] });
							}
						run_round(redo);

						passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
						++num_resolve_attempts;
					}
					midpath_report_.num_resolve_attempts = num_resolve_attempts;
					midpath_report_.passed = passed;
					if (!passed)
						std::cerr << "warning: " << midpath_report_.num_crossings_detected
						          << " path crossing(s) detected at the endgame boundary remained unresolved "
						          << "after " << num_resolve_attempts << " re-track attempt(s); "
						          << "the affected solutions may be wrong.  Consider a higher-order predictor "
						          << "or a tighter tracking tolerance.  See EndgameBoundaryMetadata()." << std::endl;

					PostEGAction();
				}
				else // worker
				{
					auto run_worker_round = [&]()
					{
						if (n_threads <= 1)
						{
							// Single-threaded worker: execute whole paths on the rank's member
							// tracker/endgame directly, from rank 0's authoritative start point.
							parallel::RunWorkerLoop<Task, Result>(comm,
								[this](Task const& task){ ExecuteOnePath(MemberDuringEGContext(), static_cast<SolnIndT>(task.path_index), task.start_point); },
								[this](Task const& task) -> Result { return PackFullPathResult(static_cast<SolnIndT>(task.path_index)); });
						}
						else
						{
							// Threaded worker: each thread owns a self-contained PathThreadState clone
							// so concurrent whole-path execution shares no mutable state.  Cloned per
							// round so it inherits any escalation (predictor bump) applied to the rank's
							// member tracker between rounds.
							auto state_factory = [this]() -> std::unique_ptr<PathThreadState>
							{
								// Clone(), not copy: System copies are SHALLOW (shared expression-tree
								// nodes whose value caches mutate on Eval); Clone() deep-copies the tree.
								// Trackers/Endgames have no default ctor, so aggregate-init via new.
								std::unique_ptr<PathThreadState> s(
									new PathThreadState{ Clone(GetTracker().GetSystem()), Clone(TargetSystem()),
									                     GetTracker(), GetEndgame(), {}, {} });
								s->tracker.SetSystem(s->sys);
								s->endgame.SetTracker(s->tracker);
								SetThreadPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
								return s;
							};

							auto track_fn = [this](std::unique_ptr<PathThreadState>& state, Task const& task) -> Result
							{
								auto idx = static_cast<SolnIndT>(task.path_index);
								ExecuteOnePath(state->Context(), idx, task.start_point);
								return PackFullPathResult(idx);
							};

							parallel::RunWorkerLoopThreaded<Task, Result>(
								comm, state_factory, track_fn, n_threads);
						}
					}; // run_worker_round

					while (true)
					{
						run_worker_round();

						int keep_going = 0;
						MPI_Bcast(&keep_going, 1, MPI_INT, 0, comm);
						if (!keep_going)
							break;

						// Match the escalation the serial RunMidpathResolution applies between
						// re-track passes: tighten midpath_retrack_tolerance_ (read by ExecuteOnePath)
						// and bump a low-order predictor on this rank's member tracker, so the next
						// round's thread clones inherit it.
						EscalateRetrackSettings();
					}
				}
			}
#endif // BERTINI2_HAVE_MPI

			// Definitions are out-of-line at the bottom of this file, after
			// output.hpp is included (avoiding a circular-include chicken-and-egg).
			void WriteMainData(std::ostream& out) const override;
			void WriteRawData(std::ostream& out)  const override;
			void ApplyParsedConfigs(std::string const& config_str) override;

			virtual ~ZeroDim() = default;
/// setup functions


			/**
			\brief Check to ensure that target system is valid for solving.
			*/
			void ConsistencyCheck() const
			{
				if (TargetSystem().HavePathVariable())
					throw std::runtime_error("unable to perform zero dim solve on target system -- has path variable, use user homotopy instead.");

				// A square zero-dim system needs one equation per dimension.  Each projective
				// (homogeneous) variable group of size k spans P^{k-1}: its k coordinates carry
				// only k-1 dimensions because scale is free, so it needs one fewer equation than
				// it has variables.  Subtract that free scale per projective group before
				// comparing -- otherwise a genuinely square multiprojective system (e.g. the
				// eigenvalue problem (A - lam I)x = 0 with x projective, lam affine) is wrongly
				// rejected.  Affine-only systems have no hom variable groups, so this is a no-op
				// for them.  (Patches, which would also enter NumTotalFunctions, are added later
				// during system preparation; this check runs on the as-supplied system.)
				if (TargetSystem().NumVariables() - TargetSystem().NumHomVariableGroups() > TargetSystem().NumTotalFunctions())
					throw std::runtime_error("unable to perform zero dim solve on target system -- underconstrained, so has no zero dimensional solutions.");

				if (!TargetSystem().IsPolynomial())
					throw std::runtime_error("unable to perform zero dim solve on target system -- system is non-polynomial, use user homotopy instead.");
			}



			void DefaultSetup()
			{
				DefaultSettingsSetup();
				DefaultSystemSetup();
				DefaultTrackerSetup();
				DefaultMidpathSetup();
			}




			void DefaultSystemSetup()
			{
				SystemManagementPolicy::SystemSetup(this->template Get<ZeroDimConf>().path_variable_name);
				num_start_points_ = StartSystem().NumStartPoints(); // populate the internal variable
			}



			/**
			Fills the tolerances and retrack settings from default values.
			*/
			void DefaultSettingsSetup()
			{
				// this code can be made generic using Boost.Hana.  
				// see https://stackoverflow.com/questions/28764085/how-to-create-an-element-for-each-type-in-a-typelist-and-add-it-to-a-vector-in-c,
				// for example
				//
				// all that would need to be done is to extract a hana::tuple_t from the 
				// typelist contained in this::Config, and then hana::for_each() over it.

				this->template Set<Tolerances>(Tolerances());
				this->template Set<PostProcessing>(PostProcessing());
				// ZeroDimConfig is not templated on the complex type, so it cannot pick its own
				// per-type ambient-precision default.  Set it here, where the tracking type is known.
				ZeroDimConf zdc;
				zdc.initial_ambient_precision = DefaultInitialAmbientPrecision<BaseComplexT>();
				this->template Set<ZeroDimConf>(zdc);
				this->template Set<AutoRetrack>(AutoRetrack());
			}

			void SetMidpathRetrackTol(NumErrorT const& rt)
			{
				midpath_retrack_tolerance_ = rt;
			}

			const auto& MidpathRetrackTol() const
			{
				return midpath_retrack_tolerance_;
			}

			/**
			\brief Set the precision at which start points are computed.

			Defaults to the initial ambient precision.  Set this to carry a higher precision forward
			(e.g. when one solve's output seeds the next).  Once set, PreSolveSetup will not override it.
			*/
			void SetStartPointPrecision(unsigned p)
			{
				start_point_precision_ = p;
				start_point_precision_set_by_user_ = true;
			}

			unsigned StartPointPrecision() const
			{
				return start_point_precision_;
			}


			/**
			call this after setting up the tolerances, etc.
			*/
			void DefaultMidpathSetup()
			{
				midpath_ = MidpathType(MidPathConfig());
			}


			void SetMidpath(MidPathConfig const& mp)
			{
				midpath_.Set(mp);
			}

			/**
			\brief Takes the default action to set up the zero dim algorithm with the default constructed tracker.

			\note should be called after the homotopy is set up, ideally.
			*/
			void DefaultTrackerSetup()
			{
				tracker_.SetSystem(Homotopy());
				tracker_.Setup(tracking::predict::DefaultPredictor(),
				              	this->template Get<Tolerances>().newton_before_endgame,
				              	this->template Get<Tolerances>().path_truncation_threshold,
								tracking::SteppingConfig(), tracking::NewtonConfig());

				tracker_.PrecisionSetup(PrecisionConfig(Homotopy()));
			}


			/**
			\brief Sets the tracker to one you supply to this function.

			Setting the tracker associates the endgame with the tracker you pass in, too.  If this is a problem, and you need this generalized so you can use a different tracker for the pre/endgame zones, please file an issue requesting this feature.

			Assumes you have done all necessary setup to it, including associating it with the homotopy for the ZeroDim algorithm.

			Again, YOU must ensure the tracker is associated with the correct homotopy
			*/
			void SetTracker(TrackerType const& new_tracker)
			{
				tracker_ = new_tracker;
				endgame_.SetTracker(tracker_);
			}

			/**
			\brief Gets the tracker

			This version gets a const reference to it.
			*/
			const TrackerType & GetTracker() const
			{
				return tracker_;
			}

			/**
			\brief Gets the tracker

			This version gets a mutable reference to it.
			*/
			TrackerType & GetTracker()
			{
				return tracker_;
			}


///  endgame specific stuff
			/**
			\brief Sets the endgame to one you supply to this function.

			Assumes you have done all necessary setup to it, including associating it with the homotopy for the ZeroDim algorithm.

			Also assumes that you have made the tracker inside the ZeroDim algorithm be the self-same tracker object as is used for the endgame you are setting here.
			*/
			void SetEndgame(EndgameType const& new_endgame)
			{
				endgame_ = new_endgame;
			}

			/**
			\brief Gets the endgame

			This version gets a `const` reference to it.
			*/
			const EndgameType & GetEndgame() const
			{
				return endgame_;
			}

			/**
			\brief Gets the endgame

			This version gets a mutable reference to it.
			*/
			EndgameType & GetEndgame()
			{
				return endgame_;
			}


			template<typename T>
			const T & GetFromEndgame() const
			{
				return endgame_.template Get<T>();
			}

			template<typename T>
			void SetToEndgame(T const& t)
			{
				endgame_.Set(t);
			}

/// the main functions



			/**
			\brief Perform the basic Zero Dim solve algorithm.

			This function iterates over all start points to the start system, tracking from each to the endgame boundary.  At the endgame boundary, path crossings are checked for.  Multiple paths which jump onto each other will not be detected, unless you are using a certified tracker, in which case this is prevented in the first place.

			Paths which have crossed are re-run, up to a certain number of times.

			Then, the points at the endgame boundary are tracked using the prescribed endgame toward the final time.

			Finally, results are post-processed.

			It is up to you to put the output somewhere.
			*/
			void Solve()
			{
				// Reproducibility: clear this thread's mpfr constant cache (pi, etc.) at the solve
				// boundary.  mpfr caches transcendental constants at the precision last requested; the
				// Cauchy endgame needs pi (roots of unity), so a high-precision solve leaves pi cached
				// at high precision.  A subsequent lower-precision solve would otherwise reuse that
				// cached value rounded down -- differing by ~1 ULP from a freshly-computed
				// low-precision pi -- making a solve's behaviour depend on what ran before it in the
				// process (non-reproducible).  Clearing here makes each solve independent of history.
				// Cheap (one constant recompute per solve); thread-local (each worker thread clears its
				// own when it starts a solve).
				mpfr_free_cache();

				solutions_user_coords_fresh_ = false;

				PreSolveChecks();

				PreSolveSetup();

				this->NotifyObservers(AlgorithmStarted<AnyZeroDim>(*this));

				// Speculative-full-path model: carry every path all the way through (pre-endgame +
				// endgame) as one unit, then detect crossings from the collected boundary points and
				// re-run any crossed path in full.  This is the same shape the distributed solve uses
				// (one worker per whole path), so serial and distributed share the per-path primitive
				// -- including computing the start point via the same ComputeStartPoint.
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto idx = static_cast<SolnIndT>(ii);
					ExecuteOnePath(MemberDuringEGContext(), idx, ComputeStartPoint(idx));
				}

				RunMidpathResolution([this](SolnIndT idx){ ExecuteOnePath(MemberDuringEGContext(), idx, ComputeStartPoint(idx)); });

				PostEGAction();

				this->NotifyObservers(AlgorithmComplete<AnyZeroDim>(*this));
			}




			/**
			\brief Get the computed solutions, in the coordinates of the user's
			original variables.

			The stored internal solutions are dehomogenized through the target
			system, once, lazily, into a cache; repeated calls return the cached
			container by const reference.  Assumes post-solve serial access, like
			the other accessors here.

			\see SolutionsInternalCoords
			*/
			const auto& SolutionsUserCoords() const
			{
				if (!solutions_user_coords_fresh_)
				{
					solutions_user_coords_.clear();
					solutions_user_coords_.reserve(solutions_post_endgame_.size());
					for (const auto& s : solutions_post_endgame_)
					{
						// A path that failed before or during the endgame leaves its endpoint
						// slot default-constructed (zero coordinates); keep index alignment with
						// FinalSolutionMetadata (output filters on it) by emitting an empty
						// placeholder rather than trying to dehomogenize an unset point.
						if (s.size() == 0)
							solutions_user_coords_.emplace_back();
						else
							solutions_user_coords_.push_back(this->TargetSystem().DehomogenizePoint(s));
					}
					solutions_user_coords_fresh_ = true;
				}
				return solutions_user_coords_;
			}

			/**
			\brief Get the computed solutions in the solver's internal coordinates:
			homogenized, and lying on the target system's patch.

			These are the coordinates for continuing work -- start points for further
			tracking re-using the target system's patch, refinement, etc.  Take a
			user-coordinates point back to this representation with
			System::HomogenizePoint on the target system.

			\see SolutionsUserCoords
			*/
			const auto& SolutionsInternalCoords() const
			{
				return solutions_post_endgame_;
			}

			/**
			\brief Select the solution points whose parallel metadata satisfies a predicate.

			The backing for the FiniteSolutions / RealSolutions / Singular / Nonsingular convenience
			accessors.  Returns points in user coordinates by default (pass user_coords=false for the
			solver's internal coordinates), mirroring SolutionsUserCoords / SolutionsInternalCoords.
			*/
			template<typename Pred>
			SolnCont<Vec<BaseComplexT>> SolutionsWhere(Pred pred, bool user_coords = true) const
			{
				auto const& sols = user_coords ? SolutionsUserCoords() : SolutionsInternalCoords();
				auto const& md   = FinalSolutionMetadata();
				SolnCont<Vec<BaseComplexT>> out;
				for (size_t i = 0; i < md.size() && i < sols.size(); ++i)
					if (pred(md[i]))
						out.push_back(sols[i]);
				return out;
			}

			/**
			\brief The finite solutions: successful endpoints the library calls FINITE (is_finite
			applies the configured endpoint_finite_threshold).  Includes singular, nonsingular, and
			real solutions alike.  \see RealSolutions, SingularSolutions, NonsingularSolutions
			*/
			SolnCont<Vec<BaseComplexT>> FiniteSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success == SuccessCode::Success && m.is_finite; }, user_coords);
			}

			/// \brief The real finite solutions (is_real applies the configured tolerance).
			SolnCont<Vec<BaseComplexT>> RealSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success == SuccessCode::Success && m.is_finite && m.is_real; }, user_coords);
			}

			/// \brief The nonsingular finite solutions (simple, well-conditioned roots).
			SolnCont<Vec<BaseComplexT>> NonsingularSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success == SuccessCode::Success && m.is_finite && !m.is_singular; }, user_coords);
			}

			/// \brief The singular finite solutions (multiple or ill-conditioned roots).
			SolnCont<Vec<BaseComplexT>> SingularSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success == SuccessCode::Success && m.is_finite && m.is_singular; }, user_coords);
			}

			/**
			\brief Get the metadat associated with the final computed solutions
			*/
			const auto& FinalSolutionMetadata() const
			{
				return solution_final_metadata_;
			}

			/**
			\brief Get the solutions as computed at the endgame boundary
			*/
			const auto& EndgameBoundarySolutions() const
			{
				return solutions_at_endgame_boundary_;
			}

			/**
			\brief Get the report from the midpath (path-crossing) check performed at the endgame
			boundary: how many crossings were detected, which paths, how many re-track attempts were
			made, and whether the check ultimately passed.

			\see MidpathCheckReport
			*/
			const MidpathCheckReport& EndgameBoundaryMetadata() const
			{
				return midpath_report_;
			}

			/**
			\brief A concise summary of the solve: how every path ended, by category, and whether any
			path was lost.  Computed on demand from the final solution metadata and the midpath report.

			\see SolveReport
			*/
			SolveReport Report() const
			{
				return SummarizeSolve(FinalSolutionMetadata(), EndgameBoundaryMetadata());
			}

		private:

			/**
			\brief Check that the solver functor is ready to go.
			*/
			void PreSolveChecks() const
			{
				if (num_start_points_ > solutions_at_endgame_boundary_.max_size())
					throw std::runtime_error("start system has more solutions than container for results.  I refuse to continue until this has been addressed.");
			}


			void PreSolveSetup()
			{
				auto num_as_size_t = static_cast<SolnIndT>(num_start_points_);

				solution_final_metadata_.resize(num_as_size_t);
				solutions_at_endgame_boundary_.resize(num_as_size_t);
				solutions_post_endgame_.resize(num_as_size_t);

				SetMidpathRetrackTol(this->template Get<Tolerances>().newton_before_endgame);

				// Default the start-point precision to the initial ambient precision.  A caller can
				// override it via SetStartPointPrecision (e.g. to carry a higher precision forward
				// from a previous solve whose output feeds this one).
				if (!start_point_precision_set_by_user_)
					start_point_precision_ = this->template Get<ZeroDimConf>().initial_ambient_precision;
			}

			/**
			Reference-bundle over the tracker/endgame/systems/observers a single path needs.
			Lets one set of per-path execution bodies (ExecuteBeforeEG / ExecuteDuringEG) serve
			both the serial flow (members) and the distributed worker (thread-owned clones), so
			there is exactly ONE per-path code path and serial/distributed cannot diverge.
			*/
			struct BeforeEGContext
			{
				TrackerType& tracker;
				tracking::FirstPrecisionRecorder<TrackerType>&  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType>& min_max_prec;
			};
			struct DuringEGContext
			{
				// const ref: the residual/dehomogenize/precision calls are const-callable, and the
				// RefToGiven policy hands back a const target system.  A mutable thread-owned clone
				// binds here too.
				SystemType const& target_sys;
				TrackerType& tracker;
				EndgameType& endgame;
				tracking::FirstPrecisionRecorder<TrackerType>&  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType>& min_max_prec;
			};

			BeforeEGContext MemberBeforeEGContext()
			{
				return BeforeEGContext{ GetTracker(), first_prec_rec_, min_max_prec_ };
			}
			DuringEGContext MemberDuringEGContext()
			{
				return DuringEGContext{ TargetSystem(), GetTracker(), GetEndgame(), first_prec_rec_, min_max_prec_ };
			}

			/**
			\brief Compute the start point for one path, authoritatively.

			The single source of a start point.  The serial flow calls this directly; in a distributed
			solve rank 0 calls it and ships the result to the worker (the worker never derives its own),
			so start points are identical across serial and distributed runs.  Computed at
			start_point_precision_ (defaults to the initial ambient precision; settable so a chained
			solve can carry precision forward).  StartPoint() mutates the shared start system's
			expression-tree value caches, so generation is serialized.
			*/
			Vec<BaseComplexT> ComputeStartPoint(SolnIndT soln_ind)
			{
				SetThreadPrecision(start_point_precision_);
				static std::mutex start_point_mutex;
				std::lock_guard<std::mutex> lock(start_point_mutex);
				return StartSystem().template StartPoint<BaseComplexT>(soln_ind);
			}

			/**
			\brief Track one path from the start time to the endgame boundary, against `ctx`.

			The single before-endgame body, shared by the serial flow and the distributed worker.
			Uses SetThreadPrecision (thread-local) throughout, so it is correct on the main thread
			and concurrently on worker threads alike.  Writes the boundary point + the pre-endgame
			portion of the metadata.
			*/
			void ExecuteBeforeEG(BeforeEGContext ctx, SolnIndT soln_ind, Vec<BaseComplexT> const& start_point)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind));

				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					ctx.tracker.AddObserver(ctx.first_prec_rec);
					ctx.tracker.AddObserver(ctx.min_max_prec);
				}

				auto& smd = solution_final_metadata_[soln_ind];
				smd.path_index    = soln_ind;
				smd.solution_index = soln_ind;

				auto initial_prec = this->template Get<ZeroDimConf>().initial_ambient_precision;
				// SetThreadPrecision: writes thread-local only, safe from concurrent threads.
				SetThreadPrecision(initial_prec);

				// Draw the condition-number probe direction ONCE here, for the whole track of this
				// point (pre-endgame AND endgame).  We have just reseeded this thread's RNG to a value
				// determined solely by the path index, so the probe is deterministic per path and
				// identical whether the path is tracked in a serial loop or on a distributed worker --
				// regardless of how the paths were split across ranks.  It is NOT refreshed for the
				// endgame (which would churn the RNG across its sample-circle sub-tracks), so the same
				// direction persists start to finish.
				ctx.tracker.RefreshConditionDirection();

				// The times are stored precision-free (mpq_rational); materialize them as the tracking
				// complex type at the current working precision.
				BaseComplexT t_start           ( BaseRealT(this->template Get<ZeroDimConf>().start_time) );
				BaseComplexT t_endgame_boundary( BaseRealT(this->template Get<ZeroDimConf>().endgame_boundary) );

				// The start point is supplied by the caller (computed once, authoritatively, via
				// ComputeStartPoint) rather than regenerated here.  In a distributed solve rank 0
				// computes it and sends it with the task, so a worker never derives its own -- that,
				// plus authoritative pi (issue #156), is what makes the start point identical across
				// serial and distributed runs.

				// Reset to the configured initial step size for this fresh path.  In the
				// speculative-full-path model a path's endgame runs immediately before the next
				// path's pre-endgame tracking on the same tracker, and the endgame leaves
				// ReinitializeInitialStepSize(false) with a tiny step; without restoring it here the
				// next path would crawl from the start time with that tiny step (effectively a hang).
				ctx.tracker.ReinitializeInitialStepSize(true);

				// Begin tracking at the intended ambient precision rather than the start point's
				// incidental precision.  Total-degree start points are generated at
				// LowestMultiplePrecision (the generator's arithmetic widens past the requested
				// digits, see issue #308), so without this the AMP tracker would start every
				// well-conditioned path in multiprecision and never drop to double.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					ctx.tracker.SetStartPrecision(initial_prec);

				Vec<BaseComplexT> result;
				auto tracking_success = ctx.tracker.TrackPath(result, t_start, t_endgame_boundary, start_point);

				solutions_at_endgame_boundary_[soln_ind] =
					EGBoundaryMetaDataT({ result, tracking_success, ctx.tracker.CurrentStepsize(), ctx.tracker.CurrentPrecision() });

				// Clear the start-precision override so it does not leak into the endgame, which
				// shares this tracker instance for its sample circles.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					ctx.tracker.SetStartPrecision(std::nullopt);

				smd.pre_endgame_success = tracking_success;

				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (ctx.first_prec_rec.DidPrecisionIncrease())
					{
						smd.precision_changed = true;
						smd.time_of_first_prec_increase = ctx.first_prec_rec.TimeOfIncrease();
					}
					else
					ctx.tracker.RemoveObserver(ctx.first_prec_rec);
					ctx.tracker.RemoveObserver(ctx.min_max_prec);
					using std::max;
					smd.max_precision_used =
						max(smd.max_precision_used, ctx.min_max_prec.MaxPrecision());
				}
			}

#ifdef BERTINI2_HAVE_MPI
			/**
			Self-contained per-thread state for executing one WHOLE path on a threaded MPI worker.
			Each std::thread owns one: a homotopy copy (tracked by `tracker`), a target-system copy
			(residual evaluation mutates System precision state), a Tracker and Endgame, and its own
			precision observers so observer-derived metadata matches serial runs.  Build a
			DuringEGContext from it to drive ExecuteOnePath.
			*/
			struct PathThreadState
			{
				System          sys;         // homotopy, tracked by `tracker`
				System          target_sys;  // for function residuals / dehomogenization
				TrackerType     tracker;
				EndgameType     endgame;
				tracking::FirstPrecisionRecorder<TrackerType>  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec;

				DuringEGContext Context()
				{
					return DuringEGContext{ target_sys, tracker, endgame, first_prec_rec, min_max_prec };
				}
			};
#endif // BERTINI2_HAVE_MPI

			/**
			\brief Apply one step of escalation to the tracker before re-tracking crossed paths.

			Two remedies are applied: (1) tighten the tracking tolerance, and (2) raise the ODE
			predictor to the default (RKF45) if a low-order predictor is in use -- a too-low-order
			predictor (notably Euler) is the most common cause of spurious boundary crossings, and
			no amount of tolerance-tightening on a first-order predictor is as effective as moving to
			a higher-order one.  The predictor bump is idempotent (once at RKF45's order it is a
			no-op).

			\note Future work: this single step wants to become a *configurable, ordered remedy list*
			(tighten tolerance, more Newton iterations, bump predictor, shrink min step size, ...),
			applied in sequence until exhausted.  That strategy abstraction is intentionally deferred;
			for now the escalation is a fixed, minimal two-remedy step.  Termination is guaranteed
			regardless, because EGBoundaryAction bounds the number of attempts.
			*/
			void EscalateRetrackSettings()
			{
				midpath_retrack_tolerance_ *= this->template Get<AutoRetrack>().midpath_decrease_tolerance_factor;
				GetTracker().SetTrackingTolerance(midpath_retrack_tolerance_);

				const auto default_predictor = tracking::predict::DefaultPredictor();
				if (tracking::predict::Order(GetTracker().GetPredictor())
				    < tracking::predict::Order(default_predictor))
					GetTracker().SetPredictor(default_predictor);
			}



			/**
			\brief Execute one whole path: start -> endgame boundary -> target, against `ctx`.

			The unit of work in the speculative-full-path model: a single worker/thread carries a
			path through both the pre-endgame tracking and the endgame, so the boundary state never
			leaves the executing context (no serialize-the-handoff, hence no precision-transfer gap).
			The pre-endgame tracking tolerance is `midpath_retrack_tolerance_` (== newton_before_endgame
			on the first pass, tightened on re-track passes by EscalateRetrackSettings); the endgame
			uses newton_during_endgame.
			*/
			void ExecuteOnePath(DuringEGContext ctx, SolnIndT soln_ind, Vec<BaseComplexT> const& start_point)
			{
				this->NotifyObservers(PathStarted<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind)));

				ctx.tracker.SetTrackingTolerance(midpath_retrack_tolerance_);
				ExecuteBeforeEG(BeforeEGContext{ ctx.tracker, ctx.first_prec_rec, ctx.min_max_prec }, soln_ind, start_point);

				if (solution_final_metadata_[soln_ind].pre_endgame_success != SuccessCode::Success)
				{
					this->NotifyObservers(PathComplete<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind)));
					return;
				}

				ctx.tracker.SetTrackingTolerance(this->template Get<Tolerances>().newton_during_endgame);
				ExecuteDuringEG(ctx, soln_ind);

				this->NotifyObservers(PathComplete<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind)));
			}


			/**
			\brief Detect path crossings from the collected boundary points and, for each crossed
			path, re-run it via `redo_one` (a full-path redo) with escalated settings.

			Shared by the serial flow and the distributed manager.  Populates midpath_report_.
			Bounded by max_num_crossed_path_resolve_attempts (0 = detect-and-report only), so it always
			terminates.  Unresolved crossings stay tracked but observable (report.passed == false +
			warning).
			*/
			template<typename RedoOne>
			void RunMidpathResolution(RedoOne&& redo_one)
			{
				auto passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());

				midpath_report_ = MidpathCheckReport{};
				midpath_report_.passed = passed;
				for (auto const& v : midpath_.GetCrossedPaths())
					midpath_report_.crossed_path_indices.push_back(v.index());
				midpath_report_.num_crossings_detected =
					static_cast<unsigned>(midpath_report_.crossed_path_indices.size());

				const auto max_attempts = this->template Get<ZeroDimConf>().max_num_crossed_path_resolve_attempts;
				unsigned num_resolve_attempts = 0;
				while (!passed && num_resolve_attempts < max_attempts)
				{
					EscalateRetrackSettings();
					for (auto const& v : midpath_.GetCrossedPaths())
						if (v.rerun())
							redo_one(static_cast<SolnIndT>(v.index()));
					passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
					++num_resolve_attempts;
				}

				midpath_report_.num_resolve_attempts = num_resolve_attempts;
				midpath_report_.passed = passed;

				if (!passed)
					std::cerr << "warning: " << midpath_report_.num_crossings_detected
					          << " path crossing(s) detected at the endgame boundary remained unresolved "
					          << "after " << num_resolve_attempts << " re-track attempt(s); "
					          << "the affected solutions may be wrong.  Consider a higher-order predictor "
					          << "or a tighter tracking tolerance.  See EndgameBoundaryMetadata()." << std::endl;
			}


			/**
			\brief Run the endgame on one path (boundary→target), against `ctx`.

			The single during-endgame body, shared by the serial flow and the distributed worker.
			Resumes at the precision the path was using at the boundary
			(`solutions_at_endgame_boundary_[idx].precision`) and writes the final solution + the
			endgame portion of the metadata.
			*/
			void ExecuteDuringEG(DuringEGContext ctx, SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind) + static_cast<uint64_t>(num_start_points_));

				auto& smd = solution_final_metadata_[soln_ind];
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (!smd.precision_changed)
						ctx.tracker.AddObserver(ctx.first_prec_rec);
					ctx.tracker.AddObserver(ctx.min_max_prec);
				}

				const auto& bdry_point = solutions_at_endgame_boundary_[soln_ind].path_point;

				ctx.tracker.SetStepSize(solutions_at_endgame_boundary_[soln_ind].last_used_stepsize);
				ctx.tracker.ReinitializeInitialStepSize(false);

				// Resume the endgame at the precision the path was actually using at the boundary,
				// rather than inferring it from the (always-multiprecision) point's mantissa, which
				// is unreliable for paths that tracked in double.
				auto start_prec = solutions_at_endgame_boundary_[soln_ind].precision;
				SetThreadPrecision(start_prec);

				ctx.endgame.SetBoundaryTime(BaseComplexT(BaseRealT(this->template Get<ZeroDimConf>().endgame_boundary)));
				ctx.endgame.SetTargetTime  (BaseComplexT(BaseRealT(this->template Get<ZeroDimConf>().target_time)));

				auto eg_success = ctx.endgame.Run(bdry_point);

				solutions_post_endgame_[soln_ind] = ctx.endgame.template FinalApproximation<BaseComplexT>();

				smd.endgame_success = eg_success;

				// an unsuccessful endgame has no final approximation, so the final-point-dependent
				// metadata cannot be computed.
				if (eg_success != SuccessCode::Success)
				{
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						ctx.tracker.RemoveObserver(ctx.first_prec_rec);
						ctx.tracker.RemoveObserver(ctx.min_max_prec);
					}
					return;
				}
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (!smd.precision_changed)
					{
						if (ctx.first_prec_rec.DidPrecisionIncrease())
						{
							smd.precision_changed = true;
							smd.time_of_first_prec_increase = ctx.first_prec_rec.TimeOfIncrease();
						}
						ctx.tracker.RemoveObserver(ctx.first_prec_rec);
					}
					ctx.tracker.RemoveObserver(ctx.min_max_prec);
					using std::max;
					smd.max_precision_used =
						max(smd.max_precision_used, ctx.min_max_prec.MaxPrecision());
				}
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					assert(Precision(solutions_post_endgame_[soln_ind])==Precision(ctx.endgame.template FinalApproximation<BaseComplexT>()));
					SetThreadPrecision(Precision(solutions_post_endgame_[soln_ind]));
					ctx.target_sys.precision(Precision(solutions_post_endgame_[soln_ind]));
				}
				smd.function_residual = static_cast<NumErrorT>(ctx.target_sys.Eval(solutions_post_endgame_[soln_ind]).template lpNorm<Eigen::Infinity>());
				smd.final_time_used = ctx.endgame.LatestTime();
				smd.condition_number = ctx.tracker.LatestConditionNumber();
				smd.newton_residual = ctx.tracker.LatestNormOfStep();

				smd.accuracy_estimate = ctx.endgame.ApproximateError();
				smd.accuracy_estimate_user_coords =
					static_cast<NumErrorT>( (ctx.target_sys.DehomogenizePoint(solutions_post_endgame_[soln_ind]) -
					ctx.target_sys.DehomogenizePoint(ctx.endgame.template PreviousApproximation<BaseComplexT>())).template lpNorm<Eigen::Infinity>() );
				smd.cycle_num = ctx.endgame.CycleNumber();
			}


			void PostEGAction()
			{
				ComputePostTrackMetadata();
			}



			/**
			\brief Populates the result_ member, with the post-processed results of the zero dim solve.

			output: result_ member variable.
			*/
			void ComputePostTrackMetadata()
			{
				ClassifyFiniteAndReal();
				ComputeMultiplicities();
				ClassifySingular();
			}

			/**
			\brief Classify each successful endpoint as finite/infinite, and (when finite) real.

			Uses the same dehomogenize-then-infinity-norm measurement the endgames use for their
			`Security::max_norm` divergence test (`System::InfinityNormOfDehomogenized`), compared
			against `endpoint_finite_threshold` -- so the metadata can never contradict the
			endgame's own verdict.  Endpoints the endgame already flagged as diverging
			(`GoingToInfinity` / `SecurityMaxNormReached`) are taken as infinite without
			recomputation.  A finite endpoint is real if the infinity norm of the imaginary parts
			of its dehomogenized coordinates is below `real_threshold`.
			*/
			void ClassifyFiniteAndReal()
			{
				using std::abs; using std::imag;
				const auto& post = this->template Get<PostProcessing>();

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto& smd = solution_final_metadata_[ii];

					if (smd.endgame_success==SuccessCode::GoingToInfinity ||
					    smd.endgame_success==SuccessCode::SecurityMaxNormReached)
					{
						smd.is_finite = false; // the endgame already decided this path diverges
						continue;
					}
					if (smd.endgame_success!=SuccessCode::Success)
						continue;              // failed otherwise: leave defaults (not finite/real/singular)

					auto user_pt = this->TargetSystem().DehomogenizePoint(solutions_post_endgame_[ii]);

					smd.is_finite =
						static_cast<NumErrorT>(user_pt.template lpNorm<Eigen::Infinity>()) <= post.endpoint_finite_threshold;

					if (smd.is_finite)
					{
						NumErrorT max_imag{0};
						for (Eigen::Index k{0}; k < user_pt.size(); ++k)
						{
							NumErrorT a = static_cast<NumErrorT>(abs(imag(user_pt(k))));
							if (a > max_imag) max_imag = a;
						}
						smd.is_real = max_imag < post.real_threshold;
					}
				}
			}

			/**
			\brief Cluster identical endpoints to compute multiplicities.

			Two endpoints are the same point when the infinity norm of the difference of their
			*dehomogenized* coordinates is below `final_tolerance * same_point_tolerance_multiplier`.
			Comparing dehomogenized (user) coordinates -- not the internal homogenized, on-patch
			coordinates, which carry the homogenizing variable and patch scaling -- is essential;
			the infinity norm matches the convergence norm the endgames use.  Only finite,
			successful endpoints are clustered (an at-infinity endpoint has no meaningful
			dehomogenized coordinates to compare).
			*/
			void ComputeMultiplicities()
			{
				const NumErrorT same_tol =
					this->template Get<Tolerances>().final_tolerance *
					this->template Get<PostProcessing>().same_point_tolerance_multiplier;

				std::vector<Vec<BaseComplexT>> user_pts(num_start_points_);
				std::vector<char> eligible(num_start_points_, 0);
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto const& smd = solution_final_metadata_[ii];
					if (smd.endgame_success==SuccessCode::Success && smd.is_finite)
					{
						user_pts[ii] = this->TargetSystem().DehomogenizePoint(solutions_post_endgame_[ii]);
						eligible[ii] = 1;
					}
				}

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					if (!eligible[ii]) continue;
					for (decltype(num_start_points_) jj{ii+1}; jj < num_start_points_; ++jj)
					{
						if (!eligible[jj]) continue;
						if ( static_cast<NumErrorT>((user_pts[ii] - user_pts[jj]).template lpNorm<Eigen::Infinity>()) < same_tol )
						{
							++solution_final_metadata_[ii].multiplicity;
							++solution_final_metadata_[jj].multiplicity;
						}
					}
				}
			}

			/**
			\brief Classify each successful endpoint as singular or not.

			Matching Bertini 1: an endpoint is singular if it is the endpoint of multiple paths
			(multiplicity > 1), or if the approximation of its condition number (spectral norm, as
			estimated by the tracker) exceeds `condition_number_threshold`.
			*/
			void ClassifySingular()
			{
				const NumErrorT cond_threshold = this->template Get<PostProcessing>().condition_number_threshold;
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto& smd = solution_final_metadata_[ii];
					if (smd.endgame_success!=SuccessCode::Success)
						continue;
					smd.is_singular = (smd.multiplicity > 1) || (smd.condition_number > cond_threshold);
				}
			}



		///////
		//	MPI pack/store helpers (only compiled when BERTINI2_HAVE_MPI is defined)
		///////

#ifdef BERTINI2_HAVE_MPI
			// Pack everything a worker computed for one whole path into the single result message.
			parallel::FullPathResult<BaseComplexT> PackFullPathResult(SolnIndT idx) const
			{
				parallel::FullPathResult<BaseComplexT> r;
				r.path_index          = idx;
				r.pre_endgame_success = solutions_at_endgame_boundary_[idx].success_code;
				r.boundary_point      = solutions_at_endgame_boundary_[idx].path_point;
				r.boundary_stepsize   = solutions_at_endgame_boundary_[idx].last_used_stepsize;
				r.boundary_precision  = solutions_at_endgame_boundary_[idx].precision;

				auto const& smd = solution_final_metadata_[idx];
				r.endgame_success   = smd.endgame_success;
				r.final_solution    = solutions_post_endgame_[idx];
				r.function_residual = smd.function_residual;
				r.condition_number  = smd.condition_number;
				r.newton_residual   = smd.newton_residual;
				r.final_time_used   = smd.final_time_used;
				r.accuracy_estimate = smd.accuracy_estimate;
				r.accuracy_estimate_user_coords = smd.accuracy_estimate_user_coords;
				r.cycle_num         = smd.cycle_num;
				r.precision_changed = smd.precision_changed;
				r.time_of_first_prec_increase = smd.time_of_first_prec_increase;
				r.max_precision_used = smd.max_precision_used;
				return r;
			}

			// Install a whole-path result on the manager.  Overwrites cleanly, so re-dispatching a
			// crossed path in a later resolve round simply replaces its earlier (crossed) result.
			void StoreFullPathResult(parallel::FullPathResult<BaseComplexT> const& r)
			{
				auto idx = static_cast<SolnIndT>(r.path_index);
				solutions_at_endgame_boundary_[idx] =
					EGBoundaryMetaDataT{ r.boundary_point, r.pre_endgame_success, r.boundary_stepsize, r.boundary_precision };
				solutions_post_endgame_[idx] = r.final_solution;

				auto& smd = solution_final_metadata_[idx];
				smd.path_index          = idx;
				smd.solution_index      = idx;
				smd.pre_endgame_success = r.pre_endgame_success;
				smd.endgame_success     = r.endgame_success;
				smd.function_residual   = r.function_residual;
				smd.condition_number    = r.condition_number;
				smd.newton_residual     = r.newton_residual;
				smd.final_time_used     = r.final_time_used;
				smd.accuracy_estimate   = r.accuracy_estimate;
				smd.accuracy_estimate_user_coords = r.accuracy_estimate_user_coords;
				smd.cycle_num           = r.cycle_num;
				smd.precision_changed   = r.precision_changed;
				smd.time_of_first_prec_increase = r.time_of_first_prec_increase;
				smd.max_precision_used  = r.max_precision_used;
			}
#endif // BERTINI2_HAVE_MPI


		///////
		//	private data members
		///////

			unsigned long long num_start_points_;
			NumErrorT midpath_retrack_tolerance_;
			MidpathCheckReport midpath_report_; ///< populated by EGBoundaryAction; exposed via EndgameBoundaryMetadata()
			unsigned start_point_precision_ = DoublePrecision(); ///< precision at which start points are computed; defaults to initial ambient precision (see PreSolveSetup)
			bool start_point_precision_set_by_user_ = false;


			/// observers used during tracking
			// i feel like these should be factored out into some policy class which prescribes how they are used, so that the actions taken are customizable.
			tracking::FirstPrecisionRecorder<TrackerType> first_prec_rec_;
			tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec_;


			/// function objects used during the algorithm
			TrackerType tracker_;
			EndgameType endgame_;
			MidpathType midpath_;



			/// computed data
			SolnCont< EGBoundaryMetaDataT > solutions_at_endgame_boundary_; // the BaseRealT is the last used stepsize
			SolnCont<Vec<BaseComplexT> > solutions_post_endgame_;
			mutable SolnCont<Vec<BaseComplexT> > solutions_user_coords_; ///< lazy cache for SolutionsUserCoords; serial post-solve access assumed
			mutable bool solutions_user_coords_fresh_ = false;
			SolnCont<SolutionMetaDataT> solution_final_metadata_;


		}; // struct ZeroDim

	} // ns algo

} // ns bertini

// Include output formatters after ZeroDim is fully defined.
// output.hpp includes zero_dim_solve.hpp, so #pragma once prevents re-inclusion
// and the circular dependency is resolved.  Any TU that gets zero_dim_solve.hpp
// therefore also gets output.hpp, making WriteMainData/WriteRawData instantiable.
#include "bertini2/nag_algorithms/output.hpp"

namespace bertini {
namespace algorithm {

template<typename TrackerType, typename EndgameType,
         typename SystemType, typename StartSystemType,
         template<typename,typename> class SystemManagementP>
inline void
ZeroDim<TrackerType,EndgameType,SystemType,StartSystemType,SystemManagementP>::WriteMainData(std::ostream& out) const
{
	output::Classic<ZeroDim>::MainData(out, *this);
}

template<typename TrackerType, typename EndgameType,
         typename SystemType, typename StartSystemType,
         template<typename,typename> class SystemManagementP>
inline void
ZeroDim<TrackerType,EndgameType,SystemType,StartSystemType,SystemManagementP>::WriteRawData(std::ostream& out) const
{
	output::Classic<ZeroDim>::RawData(out, *this);
}

} // ns algorithm
} // ns bertini

// Include config parsers after ZeroDim is fully defined, then provide the
// out-of-line definition of ApplyParsedConfigs.  Parsing headers do not include
// zero_dim_solve.hpp, so there is no circular dependency here.
#include "bertini2/io/parsing/settings_parsers.hpp"

namespace bertini {
namespace algorithm {

// Injects every element of a std::tuple<Ts...> into `target` via Set<T>.
// Works for any Configured<>-derived target whose typelist contains all Ts.
template<typename Target, typename... Ts>
void InjectParsedTuple(Target& target, std::tuple<Ts...> const& t) {
	(target.template Set<Ts>(std::get<Ts>(t)), ...);
}

template<typename TrackerType, typename EndgameType,
         typename SystemType, typename StartSystemType,
         template<typename,typename> class SystemManagementP>
void
ZeroDim<TrackerType,EndgameType,SystemType,StartSystemType,SystemManagementP>
    ::ApplyParsedConfigs(std::string const& config_str)
{
	using namespace parsing::classic;

	// 1. ZeroDim-owned configs (Tolerances, PostProcessing, ZeroDimConf, AutoRetrack)
	using ZDConfs = typename Config::UsedConfigs;
	auto zd = ConfigParser<ZDConfs>::Parse(config_str);
	InjectParsedTuple(*this, zd);
	// The path variable name comes from ZeroDimConf, so the homotopy must be
	// re-formed if the parsed name differs from the one used at construction.
	// Only re-run setup in that case: re-running unconditionally re-randomizes
	// gamma and the start system for no reason, and (before Homogenize() was made
	// idempotent) re-corrupted the prepared target system.
	if (!Homotopy().HavePathVariable()
	    || Homotopy().GetPathVariable()->name() != this->template Get<ZeroDimConf>().path_variable_name)
	{
		DefaultSystemSetup();
	}

	// 2. Tracker — uses positional Setup() rather than Set<T>, so handle explicitly.
	using TkConfs = detail::TypeList<
	    tracking::SteppingConfig,
	    tracking::NewtonConfig,
	    tracking::Predictor>;
	auto tk = ConfigParser<TkConfs>::Parse(config_str);
	tracker_.Setup(
	    std::get<tracking::Predictor>(tk),
	    this->template Get<Tolerances>().newton_before_endgame,
	    this->template Get<Tolerances>().path_truncation_threshold,
	    std::get<tracking::SteppingConfig>(tk),
	    std::get<tracking::NewtonConfig>(tk));
	tracker_.PrecisionSetup(PrecisionConfig(Homotopy()));
	endgame_.SetTracker(tracker_); // keep endgame's tracker ref consistent after reconfigure

	// 3. Endgame-owned configs — endgame_ is the concrete EndgameType deriving from
	//    Configured<its_configs>, so Set<T> is available directly.
	//    AlgoTraits<EndgameType>::NeededConfigs is public (unlike EndgameType::Configs).
	using EGConfs = typename endgame::AlgoTraits<EndgameType>::NeededConfigs;
	auto eg = ConfigParser<EGConfs>::Parse(config_str);
	InjectParsedTuple(endgame_, eg);

	// 4. Midpath config
	SetMidpath(ConfigParser<MidPathConfig>::Parse(config_str));
}

} // ns algorithm
} // ns bertini


// Explicit instantiation declarations — suppress re-instantiation of the six
// production ZeroDim types in every including TU.  Definitions live in
// core/src/eti/zero_dim_eti.cpp; see ADR-0014.  Other combos (different start
// systems, RefToGiven policy) simply instantiate implicitly as before.
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"

namespace bertini{ namespace algorithm{

extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::PSEG,     System, start_system::TotalDegree>;
extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::Cauchy,   System, start_system::TotalDegree>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::PSEG,   System, start_system::TotalDegree>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::Cauchy, System, start_system::TotalDegree>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::PSEG,                 System, start_system::TotalDegree>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::Cauchy,               System, start_system::TotalDegree>;

// the blackbox switch ladder additionally reaches MHomogeneous (CloneGiven) and
// User (RefToGiven) starts; definitions in core/src/eti/zero_dim_blackbox_eti.cpp
extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::PSEG,     System, start_system::MHomogeneous>;
extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::Cauchy,   System, start_system::MHomogeneous>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::PSEG,   System, start_system::MHomogeneous>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::Cauchy, System, start_system::MHomogeneous>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::PSEG,                 System, start_system::MHomogeneous>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::Cauchy,               System, start_system::MHomogeneous>;

extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::PSEG,     System, start_system::User, policy::RefToGiven>;
extern template struct ZeroDim<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::Cauchy,   System, start_system::User, policy::RefToGiven>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::PSEG,   System, start_system::User, policy::RefToGiven>;
extern template struct ZeroDim<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::Cauchy, System, start_system::User, policy::RefToGiven>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::PSEG,                 System, start_system::User, policy::RefToGiven>;
extern template struct ZeroDim<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::Cauchy,               System, start_system::User, policy::RefToGiven>;

}} // namespaces
