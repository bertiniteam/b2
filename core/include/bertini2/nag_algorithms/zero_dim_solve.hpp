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



		template<typename TrackerType, typename EndgameType, typename SystemType, typename StartSystemType>
		struct ZeroDim : public Observable<>
		{
			BERTINI_DEFAULT_VISITABLE();

			using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    = typename tracking::TrackerTraits<TrackerType>::BaseRealType;
			
			using PrecisionConfig = typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

			using SolnIndT = typename SolnCont<BaseComplexType>::size_type;
			
			using MidpathT = Midpath<StartSystemType, BaseRealType, BaseComplexType>;

			struct MetaData{

				// only vaguely metadata.  artifacts of randomness or ordering
				unsigned long long path_index;     		// path number of the solution
				unsigned long long solution_index;      	// solution number

				// computed across all of the solve
				bool precision_changed = false;
				BaseComplexType time_of_first_prec_increase;    // time value of the first increase in precision

				



				//computed in pre-endgame only

				tracking::SuccessCode pre_endgame_success;     // success code 



				// computed in endgame only
				BaseRealType condition_number; 		// the latest estimate on the condition number
				BaseRealType newton_residual; 		// the latest newton residual 

				BaseComplexType final_t;   					// the final value of time
				BaseRealType accuracy_estimate; 			// accuracy estimate between extrapolations

				unsigned cycle_num;    	// cycle number used in extrapolations
				
				tracking::SuccessCode endgame_success;      		// success code 


				// things added by post-processing
				BaseRealType function_residual; 	// the latest function residual

				int multiplicity; 		// multiplicity
				bool is_real;       		// real flag:  0 - not real, 1 - real
				bool is_finite;     		// finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
				bool is_sing;       		// singular flag: 0 - non-sigular, 1 - singular
			};


			ZeroDim(System const& sys) : target_system_(Clone(sys)), tracker_(target_system_), endgame_(tracker_)
			{
				ConsistencyCheck();
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
				using Variable = node::Variable;
				using Var = std::shared_ptr<Variable>;

				target_system_.Homogenize(); // work over projective coordinates
				target_system_.AutoPatch(); // then patch if needed

				start_system_ = StartSystemType(target_system_); // make the start system from the target system.

				num_start_points_ = start_system_.NumStartPoints();

				Var t = std::make_shared<Variable>("ZERO_DIM_PATH_VARIABLE"); 

				homotopy_ = (1-t)*target_system_ + std::make_shared<node::Rational>(node::Rational::Rand())*t*start_system_;
				homotopy_.AddPathVariable(t);
				assert(homotopy_.IsHomogeneous());

				tolerances_ = algorithm::config::Tolerances<BaseRealType>();
				
				boundary_near_tol_ = NumTraits<BaseRealType>::FromString("1e-5");
				midpath_ = std::make_shared<MidpathT>(start_system_, boundary_near_tol_);

				tracker_ = TrackerType(homotopy_);

				tracker_.Setup(tracking::predict::DefaultPredictor(),
				              	tolerances_.newton_before_endgame, NumTraits<BaseRealType>::FromString("1e5"),
								tracking::config::Stepping<BaseRealType>(), tracking::config::Newton());

				tracker_.PrecisionSetup(PrecisionConfig(homotopy_));
				
				initial_ambient_precision_ = DoublePrecision();

				t_start_ = static_cast<BaseComplexType>(1);
				t_endgame_boundary_ = NumTraits<BaseComplexType>::FromString("0.1");
				t_end_ = static_cast<BaseComplexType>(0);	

				setup_complete_ = true;
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

			bool SetupComplete() const
			{
				return setup_complete_;
			}


		private:

			/**
			\brief Check that the solver functor is ready to go.
			*/
			void PreSolveChecks()
			{
				if (!SetupComplete())
					throw std::runtime_error("attempting to Solve ZeroDim, but setup was not completed");

				if (num_start_points_ > solutions_at_endgame_boundary_.max_size())
					throw std::runtime_error("start system has more solutions than container for results.  I refuse to continue until this has been addressed.");

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
					DefaultPrecision(initial_ambient_precision_);

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

			void EGBoundaryAction()
			{
				auto midcheck = midpath_->Check(solutions_at_endgame_boundary_);
				
				unsigned num_resolve_attempts = 0;
				while (!midcheck.Passed() && num_resolve_attempts < zero_dim_config_.max_num_crossed_path_resolve_attempts)
				{
					MidpathResolve(midcheck, solutions_at_endgame_boundary_);
					midcheck = midpath_->Check(solutions_at_endgame_boundary_);
					num_resolve_attempts++;
				}
			}

			template<typename T>
			void MidpathResolve(typename MidpathT::Data const&,
			                    T const& solutions_at_endgame_boundary)
			{
				// insert code here to retrack the crossed paths.
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
			config::ZeroDim zero_dim_config_;
			config::Tolerances<BaseRealType> tolerances_;
			BaseRealType boundary_near_tol_;

			PrecisionConfig precision_config_;
			
			// do not permute the order of these System declarations
			System target_system_;
			StartSystemType start_system_;
			System homotopy_;

			TrackerType tracker_;
			EndgameType endgame_;

			BaseComplexType t_start_, t_endgame_boundary_, t_end_;

			SolnCont< std::tuple<Vec<BaseComplexType>, tracking::SuccessCode, BaseRealType>> solutions_at_endgame_boundary_;
			
			std::shared_ptr<MidpathT> midpath_;

			SolnCont<Vec<BaseComplexType> > endgame_solutions_;
			SolnCont<MetaData> solution_metadata_;

			unsigned long long num_start_points_;


			bool setup_complete_ = false;
			unsigned initial_ambient_precision_ = DoublePrecision();

		}; // struct ZeroDim


	} // algo

} // bertini



