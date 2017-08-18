//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve/algorithm.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/zero_dim_solve/algorithm.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/zero_dim_solve/algorithm.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/zero_dim_solve/algorithm.hpp

\brief Provides the algorithm for computing all zero-dimensional solutions for an algberaic system.
*/

#pragma once

#include "bertini2/num_traits.hpp"

#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"
#include "bertini2/nag_algorithms/midpath_check.hpp"
#include "bertini2/io/generators.hpp"

#include "bertini2/detail/configured.hpp"

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
	using BaseRealType = typename tracking::TrackerTraits<TrackerType>::BaseRealType;
	using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;

	using NeededConfigs = detail::TypeList<
								TolerancesConfig,
								PostProcessingConfig,
								ZeroDimConfig<BaseComplexType>,
								AutoRetrackConfig
								>;
};


struct AnyZeroDim : public virtual AnyAlgorithm
{
	virtual ~AnyZeroDim() = default;
};


/**
\brief the basic zero dim algorithm, which solves a system.
*/
		template<	typename TrackerType, typename EndgameType,
					typename SystemType, typename StartSystemType,
					template<typename,typename> class SystemManagementP>
		struct ZeroDim :
							public virtual AnyZeroDim,
							public Observable<>,
							public SystemManagementP<SystemType, StartSystemType>,
							public detail::Configured<
								typename AlgoTraits< ZeroDim<TrackerType, EndgameType, SystemType, StartSystemType, SystemManagementP>>::NeededConfigs>
		{
			BERTINI_DEFAULT_VISITABLE();

/// a bunch of using statements to reduce typing.
			using BaseComplexType 	= typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    	= typename tracking::TrackerTraits<TrackerType>::BaseRealType;

			using PrecisionConfig 	= typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

			using SolnIndT 			= typename SolnCont<BaseComplexType>::size_type;



			using SystemManagementPolicy = SystemManagementP<SystemType, StartSystemType>;

			using StoredSystemT = typename SystemManagementPolicy::StoredSystemT;
			using StoredStartSystemT = typename SystemManagementPolicy::StoredStartSystemT;


			using Config = detail::Configured<
								typename AlgoTraits<ZeroDim<TrackerType, EndgameType, SystemType, StartSystemType, SystemManagementP>>::NeededConfigs>;
			using Config::Get;


			using Tolerances = TolerancesConfig;
			using PostProcessing = PostProcessingConfig;
			using ZeroDimConf = ZeroDimConfig<BaseComplexType>;
			using AutoRetrack = AutoRetrackConfig;

/// metadata structs

			struct AlgorithmMetaData
			{
				SolnIndT number_path_failures = 0;
				SolnIndT number_path_successes = 0;
				SolnIndT number_paths_tracked = 0;

				std::chrono::system_clock::time_point start_time;
				std::chrono::microseconds elapsed_time;
			};



			struct SolutionMetaData
			{

				// only vaguely metadata.  artifacts of randomness or ordering
				SolnIndT path_index;     		// path number of the solution
				SolnIndT solution_index;      	// solution number

				///// things computed across all of the solve
				bool precision_changed = false;
				BaseComplexType time_of_first_prec_increase;    // time value of the first increase in precision
				decltype(DefaultPrecision()) max_precision_used = 0;



				///// things computed in pre-endgame only
				SuccessCode pre_endgame_success = SuccessCode::NeverStarted;     // success code





				///// things computed in endgame only
				NumErrorT condition_number; 				// the latest estimate on the condition number
				NumErrorT newton_residual; 				// the latest newton residual
				BaseComplexType final_time_used;   			// the final value of time tracked to
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
			};


			struct EGBoundaryMetaData
			{
				Vec<BaseComplexType> path_point;
				SuccessCode success_code = SuccessCode::NeverStarted;
				BaseRealType last_used_stepsize;

				EGBoundaryMetaData() = default;
				EGBoundaryMetaData(EGBoundaryMetaData const&) = default;
				EGBoundaryMetaData(Vec<BaseComplexType> const& pt, SuccessCode const& code, BaseRealType const& ss) :
					path_point(pt), success_code(code), last_used_stepsize(ss)
				{}
					
			};



// a few more using statements

			using MidpathType = MidpathChecker<BaseRealType, BaseComplexType, EGBoundaryMetaData>;

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
			{}

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
						GetTracker().AddObserver(&first_prec_rec_);
						GetTracker().AddObserver(&min_max_prec_);
					}

					auto& smd = solution_final_metadata_[soln_ind];

					smd.path_index = soln_ind;
					smd.solution_index = soln_ind;

				DefaultPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
				auto t_start = this->template Get<ZeroDimConf>().start_time;
				auto t_endgame_boundary = this->template Get<ZeroDimConf>().endgame_boundary;
				auto start_point = StartSystem().template StartPoint<BaseComplexType>(soln_ind);

				Vec<BaseComplexType> result;
				auto tracking_success = GetTracker().TrackPath(result, t_start, t_endgame_boundary, start_point);

				solutions_at_endgame_boundary_[soln_ind] = EGBoundaryMetaData({ result, tracking_success, GetTracker().CurrentStepsize() });

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
						GetTracker().RemoveObserver(&first_prec_rec_);
						GetTracker().RemoveObserver(&min_max_prec_);
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
							GetTracker().AddObserver(&first_prec_rec_);
						GetTracker().AddObserver(&min_max_prec_);
					}

				const auto& bdry_point = solutions_at_endgame_boundary_[soln_ind].path_point;


				GetTracker().SetStepSize(solutions_at_endgame_boundary_[soln_ind].last_used_stepsize);
				GetTracker().ReinitializeInitialStepSize(false);

				DefaultPrecision(Precision(bdry_point));
				// we make these fresh so they are in the correct precision to start.
				BaseComplexType t_end = this->template Get<ZeroDimConf>().target_time;
				BaseComplexType t_endgame_boundary = this->template Get<ZeroDimConf>().endgame_boundary;

				auto eg_success = GetEndgame().Run(t_endgame_boundary, bdry_point, t_end);

				solutions_post_endgame_[soln_ind] = GetEndgame().template FinalApproximation<BaseComplexType>();


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
						GetTracker().RemoveObserver(&first_prec_rec_);
						GetTracker().RemoveObserver(&min_max_prec_);
						using std::max;
						smd.max_precision_used =
							max(smd.max_precision_used, min_max_prec_.MaxPrecision());
					}
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						assert(Precision(solutions_post_endgame_[soln_ind])==Precision(GetEndgame().template FinalApproximation<BaseComplexType>()));
						DefaultPrecision(Precision(solutions_post_endgame_[soln_ind]));
						TargetSystem().precision(Precision(solutions_post_endgame_[soln_ind]));
					}
					smd.function_residual = static_cast<NumErrorT>(TargetSystem().Eval(solutions_post_endgame_[soln_ind]).template lpNorm<Eigen::Infinity>());
					smd.final_time_used = GetEndgame().LatestTime();
					smd.condition_number = GetTracker().LatestConditionNumber();
					smd.newton_residual = GetTracker().LatestNormOfStep();

					smd.accuracy_estimate = GetEndgame().template ApproximateError();
					smd.accuracy_estimate_user_coords =
						static_cast<NumErrorT>( (TargetSystem().DehomogenizePoint(solutions_post_endgame_[soln_ind]) -
						TargetSystem().DehomogenizePoint(GetEndgame().template PreviousApproximation<BaseComplexType>())).template lpNorm<Eigen::Infinity>() );
					smd.cycle_num = GetEndgame().template CycleNumber();
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
			SolnCont< EGBoundaryMetaData > solutions_at_endgame_boundary_; // the BaseRealType is the last used stepsize
			SolnCont<Vec<BaseComplexType> > solutions_post_endgame_;
			SolnCont<SolutionMetaData> solution_final_metadata_;


		}; // struct ZeroDim

	} // ns algo

} // ns bertini
