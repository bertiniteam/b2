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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/zero_dim_solve.hpp 

\brief Provides the algorithm for computing all zero-dimensional solutions for an algberaic system.  
*/


#pragma once

#include "bertini2/num_traits.hpp"
#include "bertini2/nag_algorithms/config.hpp"
#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"
#include "bertini2/nag_algorithms/midpath_check.hpp"
#include "bertini2/io/generators.hpp"


namespace bertini {

	namespace algorithm {



		/**
		This system management policy makes it so that the zero dim algorithm makes a clone of the supplied system when the algorithm is created.  The zerodim algorithm will homogenize the system (that's why you want a clone, so it leaves your original system untouched), form a start system (of your type, inferred from the template parameter for the zerodim alg), and couple the two together into the homotopy used to track.

		If you don't want it to take copies, or homogenize, etc, use a different policy.
	
		Throws if you try to set your own start system or homotopy.

		\see RefToGiven
		*/
		template<typename SystemType, typename StartSystemType>
		struct CloneGiven{

			using StoredSystemT = SystemType;
			using StoredStartSystemT = StartSystemType;

			static
			StoredSystemT AtConstruct(SystemType const& sys)
			{
				return Clone(sys);
			}


			template<typename T>
			static
			T AtSet(T const& sys)
			{
				throw std::runtime_error("when using this policy, you cannot provide your own start system or homotopy.  Either change policies, or use default behaviour");
			}


			/**
			Homogenize and patch the target system.
			*/
			static 
			void PrepareTarget(SystemType & target)
			{
				// target system came from the constructor
				target.Homogenize(); // work over projective coordinates
				target.AutoPatch(); // then patch if needed
			}

			static void FormStart(StartSystemType & start, SystemType const& target)
			{
				start = StartSystemType(target);	
			}

			static
			void FormHomotopy(SystemType & homotopy, SystemType const& start, StartSystemType const& target, std::string const& path_variable_name)
			{
				auto t = MakeVariable(path_variable_name); 

				homotopy = (1-t)*target + MakeRational(node::Rational::Rand())*t*start;
				homotopy.AddPathVariable(t);
			}
		};


		/**
		This system management policy allows the user to prevent the zero dim algorithm from making clones, and instead the burden of supplying the target system, start system, and homotopy are entirely up to the user.

		Using this policy implies the user manages these things entirely.

		\see CloneGiven
		*/
		template<typename SystemType, typename StartSystemType>
		struct RefToGiven{

			using StoredSystemT = std::reference_wrapper<SystemType>;
			using StoredStartSystemT = std::reference_wrapper<StartSystemType>;

			static
			StoredSystemT AtConstruct(SystemType const& sys)
			{
				return std::ref(sys);
			}

			template<typename T>
			static
			T AtSet(T const& sys)
			{
				return std::ref(sys);
			}

			/**
			 don't do anything to the target system.  it's assumed it was passed in ready to track.
			*/
			static
			void PrepareTarget(SystemType const& sys)
			{ 
			}

			static void FormStart(StoredSystemT & start, SystemType const& target)
			{
				start = StartSystemType(target);	
			}

			static
			void FormHomotopy(SystemType & homotopy, SystemType const& target, StartSystemType const& start, std::string const& path_variable_name)
			{
				auto t = MakeVariable(path_variable_name); 

				homotopy = (1-t)*target + MakeRational(node::Rational::Rand())*t*start;
				homotopy.AddPathVariable(t);
			}

		};




		template<	typename TrackerType, typename EndgameType, 
					typename SystemType, typename StartSystemType, 
					typename SystemManagementPolicy = CloneGiven<SystemType, StartSystemType> >
		struct ZeroDim : public Observable<>
		{
			BERTINI_DEFAULT_VISITABLE();

			using BaseComplexType 	= typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    	= typename tracking::TrackerTraits<TrackerType>::BaseRealType;
			
			using PrecisionConfig 	= typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

			using SolnIndT 			= typename SolnCont<BaseComplexType>::size_type;
			
			using MidpathT 			= Midpath<StartSystemType, BaseRealType, BaseComplexType>;

			using StoredSystemT = typename SystemManagementPolicy::StoredSystemT;
			using StoredStartSystemT = typename SystemManagementPolicy::StoredStartSystemT;

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




				///// things computed in pre-endgame only
				tracking::SuccessCode pre_endgame_success;     // success code 





				///// things computed in endgame only
				BaseRealType condition_number; 				// the latest estimate on the condition number
				BaseRealType newton_residual; 				// the latest newton residual 
				BaseComplexType final_time_used;   			// the final value of time tracked to
				BaseRealType accuracy_estimate; 			// accuracy estimate between extrapolations
				unsigned cycle_num;    						// cycle number used in extrapolations
				tracking::SuccessCode endgame_success;      // success code 




				///// things added by post-processing
				BaseRealType function_residual; 	// the latest function residual

				int multiplicity; 		// multiplicity
				bool is_real;       		// real flag:  0 - not real, 1 - real
				bool is_finite;     		// finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
				bool is_singular;       		// singular flag: 0 - non-sigular, 1 - singular
			};


			/**
			Important note: pay attention to the system management policy.
			*/
			ZeroDim(System const& sys) : target_system_(SystemManagementPolicy::AtConstruct(sys)), tracker_(target_system_), endgame_(tracker_)
			{
				ConsistencyCheck();
			}


			/**
			A getter for the system to be tracked to.
			*/
			const System& TargetSystem() const
			{
				return target_system_;
			}

			/**
			A getter for the homotopy being used.
			*/
			const System& Homotopy() const
			{
				return homotopy_;
			}

			/**
			A getter for the start system being used.
			*/
			const StartSystemType & StartSystem() const
			{
				return start_system_;
			}



			/**
			\brief Check to ensure that target system is valid for solving.
			*/
			void ConsistencyCheck() const
			{
				if (target_system_.HavePathVariable())
					throw std::runtime_error("unable to perform zero dim solve on target system -- has path variable, use user homotopy instead.");

				if (target_system_.NumVariables() > target_system_.NumTotalFunctions())
					throw std::runtime_error("unable to perform zero dim solve on target system -- underconstrained, so has no zero dimensional solutions.");

				if (!target_system_.IsPolynomial())
					throw std::runtime_error("unable to perform zero dim solve on target system -- system is non-polynomial, use user homotopy instead.");
			}



			void DefaultSetup()
			{
				DefaultSettingsSetup();
				
				DefaultTimeSetup();

				DefaultSystemSetup();
				DefaultTrackerSetup();
				DefaultMidpathSetup();
				
				setup_complete_ = true;
			}


			
			/**
			Fills the tolerances and retrack settings from default values.
			*/
			void DefaultSettingsSetup()
			{
				zero_dim_config_ = algorithm::config::ZeroDim<BaseComplexType>();
				tolerances_ = algorithm::config::Tolerances<BaseRealType>();
				retrack_ = algorithm::config::AutoRetrack<BaseRealType>();
				post_processing_ = algorithm::config::PostProcessing<BaseRealType>();

			}


			/**
			call this after setting up the tolerances, etc.
			*/
			void DefaultMidpathSetup()
			{
				midpath_reduced_tolerance_ = tolerances_.newton_before_endgame;
				midpath_ = std::make_shared<MidpathT>(start_system_, retrack_.boundary_near_tol);
			}


			/**
			gets the start, target, and endgame boundary times from the already-set zero dim config.  

			please ensure the zero dim config is already set.
			*/
			void DefaultTimeSetup()
			{
				t_start_ 			= zero_dim_config_.start_time;
				t_endgame_boundary_ = zero_dim_config_.endgame_boundary;
				t_end_ 				= zero_dim_config_.target_time;
			}

			


			/**
			\brief Sets up the homotopy for the system to be solved.
			
			1. Homogenizes the system, 
			2. patches it, 
			3. constructs the start system,
			4. stores the number of start points, 
			5. makes a path variable,
			6. forms the straight line homotopy between target and start, with the gamma trick

			boom, you're ready to go.
			*/
			void DefaultSystemSetup()
			{
				SystemManagementPolicy::PrepareTarget(target_system_);

				// now we populate the start system
				SystemManagementPolicy::FormStart(start_system_, target_system_);

				num_start_points_ = start_system_.NumStartPoints(); // populate the internal variable

				SystemManagementPolicy::FormHomotopy(homotopy_, start_system_, target_system_, zero_dim_config_.path_variable_name);
			}


			/**
			\brief Set the start system and homotopy.

			Uses the management policy you put in as a template parameter when making the ZeroDim object.  Pay attention.
			*/
			void SetStartAndHomotopy(StartSystemType const& ss, System const& hh)
			{
				start_system_ = SystemManagementPolicy::AtSet(ss);
				homotopy_ = SystemManagementPolicy::AtSet(hh);
			}

			/**
			\brief Takes the default action to set up the zero dim algorithm with the default constructed tracker.

			\note should be called after the homotopy is set up, ideally.
			*/
			void DefaultTrackerSetup()
			{
				tracker_ = TrackerType(homotopy_);

				tracker_.Setup(tracking::predict::DefaultPredictor(),
				              	tolerances_.newton_before_endgame, 
				              	tolerances_.path_truncation_threshold,
								tracking::config::Stepping<BaseRealType>(), tracking::config::Newton());
				
				tracker_.PrecisionSetup(PrecisionConfig(homotopy_));
			}


			/**
			\brief Sets the tracker to one you supply to this function.

			Assumes you have done all necessary setup to it, including associating it with the homotopy for the ZeroDim algorithm.
			*/
			void SetTracker(TrackerType const& new_tracker)
			{
				tracker_ = new_tracker;
			}



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
				
				PostProcess();
			}

			/**
			\brief Outputs the finalized solutions to a target destination
			*/
			template<class GeneratorT = generators::Classic>
			void Output() const
			{
				
			}

			bool IsSetupComplete() const
			{
				return setup_complete_;
			}


		private:

			/**
			\brief Check that the solver functor is ready to go.
			*/
			void PreSolveChecks()
			{
				if (!IsSetupComplete())
					throw std::runtime_error("attempting to Solve ZeroDim, but setup was not completed");

				if (num_start_points_ > solutions_at_endgame_boundary_.max_size())
					throw std::runtime_error("start system has more solutions than container for results.  I refuse to continue until this has been addressed.");
			}


			void PreSolveSetup()
			{
				auto num_as_size_t = static_cast<SolnIndT>(num_start_points_);

				solution_metadata_.resize(num_as_size_t);
				solutions_at_endgame_boundary_.resize(num_as_size_t);
				endgame_solutions_.resize(num_as_size_t);
			}

			/**
			\brief Track from the start point in time, from each start point of the start system, to the endgame boundary.  

			Results are accumulated into an internally stored variable, solutions_at_endgame_boundary_.

			The point at the endgame boundary, as well as the success flag, and the stepsize, are all stored.
			*/
			void TrackBeforeEG()
			{
				tracker_.SetTrackingTolerance(tolerances_.newton_before_endgame);

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto soln_ind = static_cast<SolnIndT>(ii);

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						tracker_.AddObserver(&first_prec_rec_);
						tracker_.AddObserver(&min_max_prec_);
					}
					DefaultPrecision(zero_dim_config_.initial_ambient_precision);

					auto start_point = start_system_.template StartPoint<BaseComplexType>(ii);

					Vec<BaseComplexType> result;
					auto tracking_success = tracker_.TrackPath(result, t_start_, t_endgame_boundary_, start_point);

					solutions_at_endgame_boundary_[soln_ind] = std::make_tuple(result, tracking_success, tracker_.CurrentStepsize());

					solution_metadata_[soln_ind].pre_endgame_success = tracking_success;

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec && first_prec_rec_.DidPrecisionIncrease())
					{
						solution_metadata_[soln_ind].precision_changed = true;
						solution_metadata_[soln_ind].time_of_first_prec_increase = first_prec_rec_.TimeOfIncrease();
					}

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						tracker_.RemoveObserver(&first_prec_rec_);
						tracker_.RemoveObserver(&min_max_prec_);
					}
				}
			}
			
			
			
			/**
			 /brief Track a single path before we reach the endgame boundary.
			*/
			void TrackSinglePathBeforeEG(SolnIndT soln_ind)
			{
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					tracker_.AddObserver(&first_prec_rec_);
					tracker_.AddObserver(&min_max_prec_);
				}
				DefaultPrecision(zero_dim_config_.initial_ambient_precision);
				
				auto start_point = start_system_.template StartPoint<BaseComplexType>(soln_ind);
				
				Vec<BaseComplexType> result;
				auto tracking_success = tracker_.TrackPath(result, t_start_, t_endgame_boundary_, start_point);
				
				solutions_at_endgame_boundary_[soln_ind] = std::make_tuple(result, tracking_success, tracker_.CurrentStepsize());
				
				solution_metadata_[soln_ind].pre_endgame_success = tracking_success;
				
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec && first_prec_rec_.DidPrecisionIncrease())
				{
					solution_metadata_[soln_ind].precision_changed = true;
					solution_metadata_[soln_ind].time_of_first_prec_increase = first_prec_rec_.TimeOfIncrease();
				}
				
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					tracker_.RemoveObserver(&first_prec_rec_);
					tracker_.RemoveObserver(&min_max_prec_);
				}
			
			}

			void EGBoundaryAction()
			{
				auto midcheckpassed = midpath_->Check(solutions_at_endgame_boundary_);
				
				unsigned num_resolve_attempts = 0;
				while (!midcheckpassed && num_resolve_attempts < zero_dim_config_.max_num_crossed_path_resolve_attempts)
				{
					MidpathResolve();
					midcheckpassed = midpath_->Check(solutions_at_endgame_boundary_);
					num_resolve_attempts++;
				}
			}

			
			void MidpathResolve()
			{
				
				midpath_reduced_tolerance_ *= retrack_.midpath_decrease_tolerance_factor;
				tracker_.SetTrackingTolerance(midpath_reduced_tolerance_);
				
				for(auto const& v : midpath_->GetCrossedPaths())
				{
					if(v.rerun())
					{
						unsigned long long index = v.index();
						auto soln_ind = static_cast<SolnIndT>(index);
						TrackSinglePathBeforeEG(soln_ind);
					}
				}

			}


			void TrackDuringEG()
			{

				tracker_.SetTrackingTolerance(tolerances_.newton_during_endgame);

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto soln_ind = static_cast<SolnIndT>(ii);

					if (solution_metadata_[soln_ind].pre_endgame_success != tracking::SuccessCode::Success)
						continue;

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec && !first_prec_rec_.DidPrecisionIncrease())
					{
						tracker_.AddObserver(&first_prec_rec_);
						tracker_.AddObserver(&min_max_prec_);
					}


					const auto& bdry_point = std::get<Vec<BaseComplexType>>(solutions_at_endgame_boundary_[soln_ind]);


					tracker_.SetStepSize(std::get<BaseRealType>(solutions_at_endgame_boundary_[soln_ind]));
					tracker_.ReinitializeInitialStepSize(false);

					DefaultPrecision(Precision(bdry_point));

					tracking::SuccessCode endgame_success = endgame_.Run(BaseComplexType(t_endgame_boundary_),bdry_point);

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec && first_prec_rec_.DidPrecisionIncrease())
					{
						solution_metadata_[soln_ind].precision_changed = true;
						solution_metadata_[soln_ind].time_of_first_prec_increase = first_prec_rec_.TimeOfIncrease();
					}

					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						tracker_.RemoveObserver(&first_prec_rec_);
						tracker_.RemoveObserver(&min_max_prec_);
					}

					solution_metadata_[soln_ind].condition_number = tracker_.LatestConditionNumber();

					endgame_solutions_[soln_ind] = endgame_.template FinalApproximation<BaseComplexType>();
				}
			}



			/**
			\brief Populates the result_ member, with the post-processed results of the zero dim solve.

			input: 

			output: result_ member variable.
			*/
			void PostProcess()
			{
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto soln_ind = static_cast<SolnIndT>(ii);
					DefaultPrecision(Precision(endgame_solutions_[soln_ind]));
					target_system_.precision(Precision(endgame_solutions_[soln_ind]));
					auto function_residuals = target_system_.Eval(endgame_solutions_[soln_ind]);
				}
			}


			tracking::FirstPrecisionRecorder<TrackerType> first_prec_rec_;
			tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec_;

			config::PostProcessing<BaseRealType> post_processing_;
			config::ZeroDim<BaseComplexType> zero_dim_config_;
			config::Tolerances<BaseRealType> tolerances_;
			config::AutoRetrack<BaseRealType> retrack_;
			BaseRealType midpath_reduced_tolerance_;

			PrecisionConfig precision_config_;
			
			// do not permute the order of these System declarations
			StoredSystemT target_system_;
			StoredStartSystemT start_system_;
			StoredSystemT homotopy_;

			TrackerType tracker_;
			EndgameType endgame_;

			BaseComplexType t_start_, t_endgame_boundary_, t_end_;

			SolnCont< std::tuple<Vec<BaseComplexType>, tracking::SuccessCode, BaseRealType>> solutions_at_endgame_boundary_;
			
			std::shared_ptr<MidpathT> midpath_;

			SolnCont<Vec<BaseComplexType> > endgame_solutions_;
			SolnCont<SolutionMetaData> solution_metadata_;

			unsigned long long num_start_points_;


			bool setup_complete_ = false;
		}; // struct ZeroDim


	} // algo

} // bertini



