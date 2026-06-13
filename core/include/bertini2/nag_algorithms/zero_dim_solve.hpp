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
#include "bertini2/nag_algorithms/common/policies.hpp"
#include "bertini2/parallel.hpp"
#include <chrono>
#include <mutex>


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
								ZeroDimConfig<BaseComplexT>,
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
	bool is_real;       		// real flag:  0 - not real, 1 - real
	bool is_finite;     		// finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
	bool is_singular;       		// singular flag: 0 - non-sigular, 1 - singular

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
			using ZeroDimConf = ZeroDimConfig<BaseComplexT>;
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
				using BeforeResult = parallel::PathBeforeEGResult<BaseComplexT>;
				using Phase2T      = parallel::Phase2Task<BaseComplexT>;
				using DuringResult = parallel::PathDuringEGResult<BaseComplexT>;

				solutions_user_coords_fresh_ = false;

				PreSolveChecks();
				PreSolveSetup();

				if (parallel::IsManager())
				{
					// ---- Phase 1: before endgame ----
					std::queue<SolnIndT> phase1_queue;
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
						phase1_queue.push(static_cast<SolnIndT>(ii));

					parallel::RunManagerLoop<SolnIndT, BeforeResult>(comm, phase1_queue,
						[this](BeforeResult const& r){ StoreBeforeEGResult(r); });

					// Midpath check + possible re-tracks run locally on rank 0.
					// Workers are idle here. This is acceptable for Phase 1 since
					// midpath crossings are rare; a future optimization could
					// redistribute re-track work.
					EGBoundaryAction();

					// ---- Phase 2: during endgame (successful paths only) ----
					std::queue<Phase2T> phase2_queue;
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
					{
						if (solution_final_metadata_[ii].pre_endgame_success == SuccessCode::Success)
						{
							auto idx = static_cast<SolnIndT>(ii);
							Phase2T task;
							task.path_index        = idx;
							task.boundary_point    = solutions_at_endgame_boundary_[idx].path_point;
							task.boundary_stepsize = solutions_at_endgame_boundary_[idx].last_used_stepsize;
							phase2_queue.push(std::move(task));
						}
					}

					parallel::RunManagerLoop<Phase2T, DuringResult>(comm, phase2_queue,
						[this](DuringResult const& r){ StoreDuringEGResult(r); });

					PostEGAction();
				}
				else
				{
					const int n_threads = parallel::WorkerThreadCount();

					// ---- Phase 1 worker ----
					if (n_threads <= 1)
					{
						parallel::RunWorkerLoop<SolnIndT, BeforeResult>(comm,
							[this](SolnIndT const& idx)
							{
								TrackSinglePathBeforeEG(idx);
							},
							[this](SolnIndT const& idx) -> BeforeResult
							{
								return PackBeforeEGResult(idx);
							});
					}
					else
					{
						// Thread-safe Phase 1: each thread owns a System copy, a Tracker
						// copy pointed at that System, and its own pair of precision
						// observers (so observer-derived metadata matches serial runs).
						// The state is heap-allocated (unique_ptr) so the System address
						// is stable after the factory returns — the tracker's
						// reference_wrapper never dangles.
						auto state_factory = [this]() -> std::unique_ptr<Phase1ThreadState>
						{
							// Clone(), not the copy constructor: System copies are SHALLOW
							// (shared_ptr expression-tree nodes), and node value caches are
							// mutated on every Eval.  Clone() deep-copies the whole tree, so
							// each thread evaluates its own.
							//
							// Trackers have no default constructor (require a System at
							// construction), so aggregate-initialize via new rather than
							// make_unique.
							std::unique_ptr<Phase1ThreadState> s(
								new Phase1ThreadState{ Clone(GetTracker().GetSystem()), GetTracker(), {}, {} });
							// Tracker copies are fully independent (value-semantic predictor,
							// corrector, and an empty observer list); just repoint the copy
							// at its own System clone.
							s->tracker.SetSystem(s->sys);
							// Same precision the serial flow sets via DefaultPrecision()
							// in TrackBeforeEG — but thread-local, from the config (the
							// global default is not reliable on worker ranks).
							SetThreadPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
							return s;
						};

						auto track_fn = [this](std::unique_ptr<Phase1ThreadState>& state, SolnIndT const& idx) -> BeforeResult
						{
							TrackSinglePathBeforeEGWith(*state, idx);
							return PackBeforeEGResult(idx);
						};

						parallel::RunWorkerLoopThreaded<SolnIndT, BeforeResult>(
							comm, state_factory, track_fn, n_threads);
					}

					// ---- Phase 2 worker (endgame) ----
					if (n_threads <= 1)
					{
						// Serial flow sets this before its endgame loop (TrackDuringEG);
						// mirror it here so worker tolerances match serial runs.
						GetTracker().SetTrackingTolerance(this->template Get<Tolerances>().newton_during_endgame);

						parallel::RunWorkerLoop<Phase2T, DuringResult>(comm,
							[this](Phase2T const& task)
							{
								auto idx = static_cast<SolnIndT>(task.path_index);
								// Install boundary data so TrackSinglePathDuringEG can read it.
								solutions_at_endgame_boundary_[idx].path_point         = task.boundary_point;
								solutions_at_endgame_boundary_[idx].last_used_stepsize = task.boundary_stepsize;
								solutions_at_endgame_boundary_[idx].success_code       = SuccessCode::Success;
								TrackSinglePathDuringEG(idx);
							},
							[this](Phase2T const& task) -> DuringResult
							{
								return PackDuringEGResult(static_cast<SolnIndT>(task.path_index));
							});
					}
					else
					{
						// Thread-safe Phase 2: each thread owns copies of the homotopy
						// (for its tracker), the target system (for residual evaluation,
						// which mutates System precision state), the tracker, the endgame
						// (rebound to the thread's tracker), and precision observers.
						auto state_factory = [this]() -> std::unique_ptr<Phase2ThreadState>
						{
							// Clone(), not copy: see the Phase 1 factory above.
							std::unique_ptr<Phase2ThreadState> s(
								new Phase2ThreadState{ Clone(GetTracker().GetSystem()), Clone(TargetSystem()),
								                       GetTracker(), GetEndgame(), {}, {} });
							s->tracker.SetSystem(s->sys);
							s->tracker.SetTrackingTolerance(this->template Get<Tolerances>().newton_during_endgame);
							s->endgame.SetTracker(s->tracker);
							SetThreadPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
							return s;
						};

						auto track_fn = [this](std::unique_ptr<Phase2ThreadState>& state, Phase2T const& task) -> DuringResult
						{
							auto idx = static_cast<SolnIndT>(task.path_index);
							// Each idx is assigned to exactly one thread — element-wise
							// writes into these pre-sized vectors race with nobody.
							solutions_at_endgame_boundary_[idx].path_point         = task.boundary_point;
							solutions_at_endgame_boundary_[idx].last_used_stepsize = task.boundary_stepsize;
							solutions_at_endgame_boundary_[idx].success_code       = SuccessCode::Success;
							TrackSinglePathDuringEGWith(*state, idx);
							return PackDuringEGResult(idx);
						};

						parallel::RunWorkerLoopThreaded<Phase2T, DuringResult>(
							comm, state_factory, track_fn, n_threads);
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
				this->template Set<ZeroDimConf>(ZeroDimConf());
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
				solutions_user_coords_fresh_ = false;

				PreSolveChecks();

				PreSolveSetup();

				TrackBeforeEG();

				EGBoundaryAction();

				TrackDuringEG();

				PostEGAction();
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
			\brief Get the metadat associated with the final computed solutions
			*/
			const auto& FinalSolutionMetadata() const
			{
				return solution_final_metadata_;
			}

			/**
			\brief Get the solutions as computed at the endgame boundary
			*/
			const auto& EndgameBoundaryData() const
			{
				return solutions_at_endgame_boundary_;
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
			}

			/**
			\brief Track from the start point in time, from each start point of the start system, to the endgame boundary.

			Results are accumulated into an internally stored variable, solutions_at_endgame_boundary_.

			The point at the endgame boundary, as well as the success flag, and the stepsize, are all stored.
			*/
			void TrackBeforeEG()
			{
				DefaultPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);

				GetTracker().SetTrackingTolerance(this->template Get<Tolerances>().newton_before_endgame);

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					TrackSinglePathBeforeEG(static_cast<SolnIndT>(ii));
				}
			}



			/**
			 /brief Track a single path before we reach the endgame boundary.
			*/
			void TrackSinglePathBeforeEG(SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind));

					// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						GetTracker().AddObserver(first_prec_rec_);
						GetTracker().AddObserver(min_max_prec_);
					}

					auto& smd = solution_final_metadata_[soln_ind];

					smd.path_index = soln_ind;
					smd.solution_index = soln_ind;

				DefaultPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
				auto t_start = this->template Get<ZeroDimConf>().start_time;
				auto t_endgame_boundary = this->template Get<ZeroDimConf>().endgame_boundary;
				auto start_point = StartSystem().template StartPoint<BaseComplexT>(soln_ind);

				// Begin tracking at the intended ambient precision rather than the start
				// point's incidental precision.  Total-degree start points are generated at
				// LowestMultiplePrecision (the generator's arithmetic widens past the requested
				// digits, see issue #308), so without this the AMP tracker would start every
				// well-conditioned path in multiprecision and never drop to double.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					GetTracker().SetStartPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);

				Vec<BaseComplexT> result;
				auto tracking_success = GetTracker().TrackPath(result, t_start, t_endgame_boundary, start_point);

				solutions_at_endgame_boundary_[soln_ind] = EGBoundaryMetaDataT({ result, tracking_success, GetTracker().CurrentStepsize(), GetTracker().CurrentPrecision() });

				// Clear the start-precision override so it does not leak into the endgame,
				// which shares this tracker instance (endgame_(tracker_)) for its sample
				// circles and must be free to begin those at the sample points' precision.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					GetTracker().SetStartPrecision(std::nullopt);
				}

					smd.pre_endgame_success = tracking_success;

					// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (first_prec_rec_.DidPrecisionIncrease())
						{
							smd.precision_changed = true;
							smd.time_of_first_prec_increase = first_prec_rec_.TimeOfIncrease();
						}
						else
						GetTracker().RemoveObserver(first_prec_rec_);
						GetTracker().RemoveObserver(min_max_prec_);
						using std::max;
						smd.max_precision_used =
							max(smd.max_precision_used, min_max_prec_.MaxPrecision());
					}


			}

#ifdef BERTINI2_HAVE_MPI
			/**
			Self-contained per-thread tracking state for Phase 1 (before-EG) work on a
			threaded MPI worker rank.  Each std::thread owns one of these: a System copy,
			a Tracker copy repointed at that System, and its own precision observers so
			observer-derived metadata is collected exactly as in serial runs.
			*/
			struct Phase1ThreadState
			{
				System          sys;
				TrackerType     tracker;
				tracking::FirstPrecisionRecorder<TrackerType>  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec;
			};

			/**
			Per-thread state for Phase 2 (endgame).  Additionally owns a target-system
			copy (residual evaluation mutates System precision state) and an Endgame copy
			rebound to the thread's tracker.
			*/
			struct Phase2ThreadState
			{
				System          sys;         // homotopy, tracked by `tracker`
				System          target_sys;  // for function residuals / dehomogenization
				TrackerType     tracker;
				EndgameT        endgame;
				tracking::FirstPrecisionRecorder<TrackerType>  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec;
			};

			/**
			\brief Track one path to the endgame boundary using thread-owned state.

			Like TrackSinglePathBeforeEG but uses the caller's tracker and observers
			(from a Phase1ThreadState) instead of the shared members.  Intended for use
			from std::thread workers.  Precision changes are thread-local only.
			*/
			void TrackSinglePathBeforeEGWith(
				Phase1ThreadState& state,
				SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind));

					// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						state.tracker.AddObserver(state.first_prec_rec);
						state.tracker.AddObserver(state.min_max_prec);
					}

				auto& smd = solution_final_metadata_[soln_ind];
				smd.path_index    = soln_ind;
				smd.solution_index = soln_ind;

				auto initial_prec = this->template Get<ZeroDimConf>().initial_ambient_precision;
				// SetThreadPrecision: writes thread-local only, safe from concurrent threads.
				SetThreadPrecision(initial_prec);

				auto t_start           = this->template Get<ZeroDimConf>().start_time;
				auto t_endgame_boundary = this->template Get<ZeroDimConf>().endgame_boundary;

				// The start system is SHARED among threads, and StartPoint() evaluates
				// expression-tree nodes (mutating their value caches), so generation is
				// serialized.  It is trivial arithmetic compared to tracking, so the
				// mutex costs nothing measurable.  Generated AFTER SetThreadPrecision
				// so the point has the same precision a serial run would give it.
				Vec<BaseComplexT> start_point;
				{
					static std::mutex start_point_mutex;
					std::lock_guard<std::mutex> lock(start_point_mutex);
					start_point = StartSystem().template StartPoint<BaseComplexT>(soln_ind);
				}

				// Begin tracking at the intended ambient precision rather than the start
				// point's incidental precision (see TrackSinglePathBeforeEG / issue #308).
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					state.tracker.SetStartPrecision(initial_prec);

				Vec<BaseComplexT> result;
				auto tracking_success = state.tracker.TrackPath(result, t_start, t_endgame_boundary, start_point);

				// Each soln_ind is unique per thread — no locking needed.
				solutions_at_endgame_boundary_[soln_ind] =
					EGBoundaryMetaDataT({ result, tracking_success, state.tracker.CurrentStepsize(), state.tracker.CurrentPrecision() });

				// Clear the start-precision override so it does not leak into the endgame
				// (which shares this tracker instance for its sample circles).
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					state.tracker.SetStartPrecision(std::nullopt);

				smd.pre_endgame_success = tracking_success;

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (state.first_prec_rec.DidPrecisionIncrease())
						{
							smd.precision_changed = true;
							smd.time_of_first_prec_increase = state.first_prec_rec.TimeOfIncrease();
						}
						else
						state.tracker.RemoveObserver(state.first_prec_rec);
						state.tracker.RemoveObserver(state.min_max_prec);
						using std::max;
						smd.max_precision_used =
							max(smd.max_precision_used, state.min_max_prec.MaxPrecision());
					}
			}

			/**
			\brief Run the endgame on one path using thread-owned state.

			Like TrackSinglePathDuringEG but uses the caller's tracker, endgame, target
			system, and observers (from a Phase2ThreadState) instead of the shared
			members.  Precision changes are thread-local only.
			*/
			void TrackSinglePathDuringEGWith(Phase2ThreadState& state, SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind) + static_cast<uint64_t>(num_start_points_));

					auto& smd = solution_final_metadata_[soln_ind];
					// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (!smd.precision_changed)
							state.tracker.AddObserver(state.first_prec_rec);
						state.tracker.AddObserver(state.min_max_prec);
					}

				const auto& bdry_point = solutions_at_endgame_boundary_[soln_ind].path_point;

				state.tracker.SetStepSize(solutions_at_endgame_boundary_[soln_ind].last_used_stepsize);
				state.tracker.ReinitializeInitialStepSize(false);

				// Resume the endgame at the precision the path was actually using at the
				// boundary, rather than inferring it from the (always-multiprecision) point's
				// mantissa, which is unreliable for paths that tracked in double.
				auto start_prec = solutions_at_endgame_boundary_[soln_ind].precision;

				SetThreadPrecision(start_prec);

				state.endgame.SetBoundaryTime(this->template Get<ZeroDimConf>().endgame_boundary);
				state.endgame.SetTargetTime  (this->template Get<ZeroDimConf>().target_time);

				auto eg_success = state.endgame.Run(bdry_point);

				solutions_post_endgame_[soln_ind] = state.endgame.template FinalApproximation<BaseComplexT>();

					// finally, store the metadata as necessary
					smd.endgame_success = eg_success;

					// an unsuccessful endgame has no final approximation, so the
					// final-point-dependent metadata cannot be computed.
					if (eg_success != SuccessCode::Success)
					{
						if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
						{
							state.tracker.RemoveObserver(state.first_prec_rec);
							state.tracker.RemoveObserver(state.min_max_prec);
						}
						return;
					}
						// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (!smd.precision_changed)
						{
							if (state.first_prec_rec.DidPrecisionIncrease())
							{
								smd.precision_changed = true;
								smd.time_of_first_prec_increase = state.first_prec_rec.TimeOfIncrease();
							}
							state.tracker.RemoveObserver(state.first_prec_rec);
						}
						state.tracker.RemoveObserver(state.min_max_prec);
						using std::max;
						smd.max_precision_used =
							max(smd.max_precision_used, state.min_max_prec.MaxPrecision());
					}
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						assert(Precision(solutions_post_endgame_[soln_ind])==Precision(state.endgame.template FinalApproximation<BaseComplexT>()));
						SetThreadPrecision(Precision(solutions_post_endgame_[soln_ind]));
						state.target_sys.precision(Precision(solutions_post_endgame_[soln_ind]));
					}
					smd.function_residual = static_cast<NumErrorT>(state.target_sys.Eval(solutions_post_endgame_[soln_ind]).template lpNorm<Eigen::Infinity>());
					smd.final_time_used = state.endgame.LatestTime();
					smd.condition_number = state.tracker.LatestConditionNumber();
					smd.newton_residual = state.tracker.LatestNormOfStep();

					smd.accuracy_estimate = state.endgame.ApproximateError();
					smd.accuracy_estimate_user_coords =
						static_cast<NumErrorT>( (state.target_sys.DehomogenizePoint(solutions_post_endgame_[soln_ind]) -
						state.target_sys.DehomogenizePoint(state.endgame.template PreviousApproximation<BaseComplexT>())).template lpNorm<Eigen::Infinity>() );
					smd.cycle_num = state.endgame.CycleNumber();
					// end metadata gathering
			}
#endif // BERTINI2_HAVE_MPI

			void EGBoundaryAction()
			{
				auto midcheckpassed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());

				unsigned num_resolve_attempts = 0;
				while (!midcheckpassed && num_resolve_attempts < this->template Get<ZeroDimConf>().max_num_crossed_path_resolve_attempts)
				{
					MidpathResolve();
					midcheckpassed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
					num_resolve_attempts++;
				}
			}


			void MidpathResolve()
			{
				ShrinkMidpathTolerance();

				for(auto const& v : midpath_.GetCrossedPaths())
				{
					if(v.rerun())
					{
						unsigned long long index = v.index();
						auto soln_ind = static_cast<SolnIndT>(index);
						TrackSinglePathBeforeEG(soln_ind);
					}
				}

			}


			void ShrinkMidpathTolerance()
			{
				midpath_retrack_tolerance_ *= this->template Get<AutoRetrack>().midpath_decrease_tolerance_factor;
				GetTracker().SetTrackingTolerance(midpath_retrack_tolerance_);
			}



			void TrackDuringEG()
			{

				GetTracker().SetTrackingTolerance(this->template Get<Tolerances>().newton_during_endgame);

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto soln_ind = static_cast<SolnIndT>(ii);

					if (solution_final_metadata_[soln_ind].pre_endgame_success != SuccessCode::Success)
						continue;

					TrackSinglePathDuringEG(soln_ind);
				}
			}


			void TrackSinglePathDuringEG(SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind) + static_cast<uint64_t>(num_start_points_));

					auto& smd = solution_final_metadata_[soln_ind];
					// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (!smd.precision_changed)
							GetTracker().AddObserver(first_prec_rec_);
						GetTracker().AddObserver(min_max_prec_);
					}

				const auto& bdry_point = solutions_at_endgame_boundary_[soln_ind].path_point;


				GetTracker().SetStepSize(solutions_at_endgame_boundary_[soln_ind].last_used_stepsize);
				GetTracker().ReinitializeInitialStepSize(false);

				// Resume the endgame at the precision the path was actually using at the
				// boundary, rather than inferring it from the (always-multiprecision) point's
				// mantissa, which is unreliable for paths that tracked in double.
				auto start_prec = solutions_at_endgame_boundary_[soln_ind].precision;

				DefaultPrecision(start_prec);

				GetEndgame().SetBoundaryTime(this->template Get<ZeroDimConf>().endgame_boundary);
				GetEndgame().SetTargetTime  (this->template Get<ZeroDimConf>().target_time);

				auto eg_success = GetEndgame().Run(bdry_point);

				solutions_post_endgame_[soln_ind] = GetEndgame().template FinalApproximation<BaseComplexT>();


					// finally, store the metadata as necessary
					smd.endgame_success = eg_success;
						// if you can think of a way to replace this `if` with something meta, please do so.
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						if (!smd.precision_changed)
						{
							if (first_prec_rec_.DidPrecisionIncrease())
							{
								smd.precision_changed = true;
								smd.time_of_first_prec_increase = first_prec_rec_.TimeOfIncrease();
							}
						}
						GetTracker().RemoveObserver(first_prec_rec_);
						GetTracker().RemoveObserver(min_max_prec_);
						using std::max;
						smd.max_precision_used =
							max(smd.max_precision_used, min_max_prec_.MaxPrecision());
					}
					// an unsuccessful endgame has no final approximation, so the
					// final-point-dependent metadata cannot be computed.
					if (eg_success != SuccessCode::Success)
						return;
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						assert(Precision(solutions_post_endgame_[soln_ind])==Precision(GetEndgame().template FinalApproximation<BaseComplexT>()));
						DefaultPrecision(Precision(solutions_post_endgame_[soln_ind]));
						TargetSystem().precision(Precision(solutions_post_endgame_[soln_ind]));
					}
					smd.function_residual = static_cast<NumErrorT>(TargetSystem().Eval(solutions_post_endgame_[soln_ind]).template lpNorm<Eigen::Infinity>());
					smd.final_time_used = GetEndgame().LatestTime();
					smd.condition_number = GetTracker().LatestConditionNumber();
					smd.newton_residual = GetTracker().LatestNormOfStep();

					smd.accuracy_estimate = GetEndgame().ApproximateError();
					smd.accuracy_estimate_user_coords =
						static_cast<NumErrorT>( (TargetSystem().DehomogenizePoint(solutions_post_endgame_[soln_ind]) -
						TargetSystem().DehomogenizePoint(GetEndgame().template PreviousApproximation<BaseComplexT>())).template lpNorm<Eigen::Infinity>() );
					smd.cycle_num = GetEndgame().CycleNumber();
					// end metadata gathering
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
				ComputeMultiplicities();
			}

			void ComputeMultiplicities()
			{
				std::vector<std::vector<int>> multiplicity_indices(num_start_points_);

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					if (solution_final_metadata_[ii].endgame_success!=SuccessCode::Success)
						continue;

					for (decltype(num_start_points_) jj{ii+1}; jj < num_start_points_; ++jj)
					{
						if (solution_final_metadata_[jj].endgame_success!=SuccessCode::Success)
							continue;

						if ( (solutions_post_endgame_[ii] - solutions_post_endgame_[jj]).norm() < this->template Get<PostProcessing>().same_point_tolerance)
						{
							multiplicity_indices[ii].push_back(static_cast<int>(jj));
							multiplicity_indices[jj].push_back(static_cast<int>(ii));
							++solution_final_metadata_[ii].multiplicity;
							++solution_final_metadata_[jj].multiplicity;
						}
					}
				}
			}



		///////
		//	MPI pack/store helpers (only compiled when BERTINI2_HAVE_MPI is defined)
		///////

#ifdef BERTINI2_HAVE_MPI
			parallel::PathBeforeEGResult<BaseComplexT> PackBeforeEGResult(SolnIndT idx) const
			{
				parallel::PathBeforeEGResult<BaseComplexT> r;
				r.path_index           = idx;
				r.pre_endgame_success  = solutions_at_endgame_boundary_[idx].success_code;
				r.boundary_point       = solutions_at_endgame_boundary_[idx].path_point;
				r.boundary_stepsize    = solutions_at_endgame_boundary_[idx].last_used_stepsize;
				r.boundary_precision   = solutions_at_endgame_boundary_[idx].precision;
				r.precision_changed    = solution_final_metadata_[idx].precision_changed;
				r.time_of_first_prec_increase = solution_final_metadata_[idx].time_of_first_prec_increase;
				r.max_precision_used   = solution_final_metadata_[idx].max_precision_used;
				return r;
			}

			void StoreBeforeEGResult(parallel::PathBeforeEGResult<BaseComplexT> const& r)
			{
				auto idx = static_cast<SolnIndT>(r.path_index);
				solutions_at_endgame_boundary_[idx] =
					EGBoundaryMetaDataT{r.boundary_point, r.pre_endgame_success, r.boundary_stepsize, r.boundary_precision};
				auto& smd = solution_final_metadata_[idx];
				smd.path_index             = idx;
				smd.solution_index         = idx;
				smd.pre_endgame_success    = r.pre_endgame_success;
				smd.precision_changed      = r.precision_changed;
				smd.time_of_first_prec_increase = r.time_of_first_prec_increase;
				smd.max_precision_used     = r.max_precision_used;
			}

			parallel::PathDuringEGResult<BaseComplexT> PackDuringEGResult(SolnIndT idx) const
			{
				parallel::PathDuringEGResult<BaseComplexT> r;
				r.path_index        = idx;
				r.final_solution    = solutions_post_endgame_[idx];
				auto const& smd     = solution_final_metadata_[idx];
				r.endgame_success   = smd.endgame_success;
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

			void StoreDuringEGResult(parallel::PathDuringEGResult<BaseComplexT> const& r)
			{
				auto idx = static_cast<SolnIndT>(r.path_index);
				solutions_post_endgame_[idx] = r.final_solution;
				auto& smd = solution_final_metadata_[idx];
				smd.endgame_success   = r.endgame_success;
				smd.function_residual = r.function_residual;
				smd.condition_number  = r.condition_number;
				smd.newton_residual   = r.newton_residual;
				smd.final_time_used   = r.final_time_used;
				smd.accuracy_estimate = r.accuracy_estimate;
				smd.accuracy_estimate_user_coords = r.accuracy_estimate_user_coords;
				smd.cycle_num         = r.cycle_num;
				// Phase 2 may have further increased precision; take the maximum.
				if (r.precision_changed && !smd.precision_changed)
				{
					smd.precision_changed = true;
					smd.time_of_first_prec_increase = r.time_of_first_prec_increase;
				}
				using std::max;
				smd.max_precision_used = max(smd.max_precision_used, r.max_precision_used);
			}
#endif // BERTINI2_HAVE_MPI


		///////
		//	private data members
		///////

			unsigned long long num_start_points_;
			NumErrorT midpath_retrack_tolerance_;


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
