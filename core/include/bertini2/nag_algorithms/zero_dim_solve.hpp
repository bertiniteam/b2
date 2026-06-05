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
#include <chrono>


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

	EGBoundaryMetaData() = default;
	EGBoundaryMetaData(EGBoundaryMetaData const&) = default;
	EGBoundaryMetaData(Vec<ComplexT> const& pt, SuccessCode const& code, RealT const& ss) :
		path_point(pt), success_code(code), last_used_stepsize(ss)
	{}
	
	bool operator==(const EGBoundaryMetaData<ComplexT> & other){
		bool result = 
			this->path_point == other.path_point
			&& this->success_code == other.success_code
			&& this->last_used_stepsize == other.last_used_stepsize
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
				Solve();
			}

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

				if (TargetSystem().NumVariables() > TargetSystem().NumTotalFunctions())
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
				PreSolveChecks();

				PreSolveSetup();

				TrackBeforeEG();

				EGBoundaryAction();

				TrackDuringEG();

				PostEGAction();
			}




			/**
			\brief Get the final computed solutions
			*/
			const auto& FinalSolutions() const
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

				auto t_start = this->template Get<ZeroDimConf>().start_time;
				auto t_endgame_boundary = this->template Get<ZeroDimConf>().endgame_boundary;

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

				Vec<BaseComplexT> result;
				auto tracking_success = GetTracker().TrackPath(result, t_start, t_endgame_boundary, start_point);

				solutions_at_endgame_boundary_[soln_ind] = EGBoundaryMetaDataT({ result, tracking_success, GetTracker().CurrentStepsize() });

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

				auto start_prec = Precision(bdry_point);

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
							multiplicity_indices[ii].push_back(jj);
							multiplicity_indices[jj].push_back(ii);
							++solution_final_metadata_[ii].multiplicity;
							++solution_final_metadata_[jj].multiplicity;
						}
					}
				}
			}



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
	DefaultSystemSetup(); // path variable name comes from ZeroDimConf; re-run after update

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
