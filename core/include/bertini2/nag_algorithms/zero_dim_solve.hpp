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

		namespace output {

			template<typename BaseComplexType>
			struct ZeroDim{
				using BaseRealType = typename Eigen::NumTraits<BaseComplexType>::Real;
				struct MetaData{
					mpz_int path_index;     		// path number of the solution
					mpz_int solution_index;      	// solution number

					BaseRealType function_residual; 	// the latest function residual
					BaseRealType condition_number; 		// the latest estimate on the condition number
					BaseRealType newton_residual; 		// the latest newton residual 

					BaseComplexType final_t;   					// the final value of time
					BaseRealType accuracy_estimate; 			// accuracy estimate between extrapolations
					BaseRealType time_of_first_prec_increase;    // time value of the first increase in precision

					unsigned cycle_num;    	// cycle number used in extrapolations
					bool success;      		// success flag 
					int multiplicity; 		// multiplicity
					bool is_real;       		// real flag:  0 - not real, 1 - real
					bool is_finite;     		// finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
					bool is_sing;       		// singular flag: 0 - non-sigular, 1 - singular
				};

				std::vector<Vec<BaseComplexType> > final_solutions_;
				std::vector<MetaData> final_metadata_;
			};
		}


		template<typename TrackerType, typename EndgameType, typename SystemType, typename StartSystemType>
		struct ZeroDim : public Observable<>
		{
			BERTINI_DEFAULT_VISITABLE();

			using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    = typename tracking::TrackerTraits<TrackerType>::BaseRealType;
			
			using PrecisionConfig = typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

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

				Var t = std::make_shared<Variable>("ZERO_DIM_PATH_VARIABLE"); 

				homotopy_ = (1-t)*target_system_ + std::make_shared<node::Rational>(node::Rational::Rand())*t*start_system_;
				homotopy_.AddPathVariable(t);

				tolerances_ = algorithm::config::Tolerances<BaseRealType>();

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

				if (start_system_.NumStartPoints() > solutions_at_endgame_boundary_.max_size())
					throw std::runtime_error("start system has more solutions than container for results.  I refuse to continue until this has been addressed.");
			}




			/**
			\brief Track from the start point in time, from each start point of the start system, to the endgame boundary.  

			Results are accumulated into an internally stored variable, solutions_at_endgame_boundary_.

			The point at the endgame boundary, as well as the success flag, and the stepsize, are all stored.
			*/
			void TrackBeforeEG()
			{
				tracker_.SetTrackingTolerance(tolerances_.newton_before_endgame);

				auto num_paths_to_track = start_system_.NumStartPoints();

				for (decltype(num_paths_to_track) ii{0}; ii < num_paths_to_track; ++ii)
				{
					DefaultPrecision(initial_ambient_precision_);

					auto start_point = start_system_.template StartPoint<BaseComplexType>(ii);

					Vec<BaseComplexType> result;
					auto tracking_success = tracker_.TrackPath(result,t_start_,t_endgame_boundary_,start_point);

					solutions_at_endgame_boundary_.push_back(std::make_tuple(result, tracking_success, tracker_.CurrentStepsize()));
				}
			}

			void EGBoundaryAction()
			{
				auto midcheck = Midpath::Check(solutions_at_endgame_boundary_);
				
				unsigned num_resolve_attempts = 0;
				while (!midcheck.Passed() && num_resolve_attempts < zero_dim_config_.max_num_crossed_path_resolve_attempts)
				{
					MidpathResolve(midcheck, solutions_at_endgame_boundary_);
					midcheck = Midpath::Check(solutions_at_endgame_boundary_);
					num_resolve_attempts++;
				}
			}

			template<typename T>
			void MidpathResolve(Midpath::Data const&, 
			                    T const& solutions_at_endgame_boundary)
			{
				// insert code here to retrack the crossed paths.
			}


			void TrackDuringEG()
			{

				tracker_.SetTrackingTolerance(tolerances_.newton_during_endgame);

				for (const auto& s : solutions_at_endgame_boundary_)
				{
					const auto& bdry_point = std::get<0>(s);
					tracker_.SetStepSize(std::get<2>(s));
					tracker_.ReinitializeInitialStepSize(false);

					DefaultPrecision(Precision(bdry_point));

					tracking::SuccessCode endgame_success = endgame_.Run(BaseComplexType(t_endgame_boundary_),bdry_point);

					endgame_solutions_.push_back(endgame_.template FinalApproximation<BaseComplexType>());
				}
			}



			/**
			\brief Populates the result_ member, with the post-processed results of the zero dim solve.

			input: 

			output: result_ member variable.
			*/
			void PostProcess()
			{
				
			}


			
			
			config::PostProcessing<BaseRealType> post_processing_;
			config::ZeroDim zero_dim_config_;
			config::Tolerances<BaseRealType> tolerances_;

			PrecisionConfig precision_config_;
			
			// do not permute the order of these System declarations
			System target_system_;
			StartSystemType start_system_;
			System homotopy_;

			TrackerType tracker_;
			EndgameType endgame_;

			BaseComplexType t_start_, t_endgame_boundary_, t_end_;

			std::vector< std::tuple<Vec<BaseComplexType>, tracking::SuccessCode, BaseRealType>> solutions_at_endgame_boundary_;
			std::vector<Vec<BaseComplexType> > endgame_solutions_;


			bool setup_complete_ = false;
			unsigned initial_ambient_precision_ = DoublePrecision();


			output::ZeroDim<BaseComplexType> results_;
		}; // struct ZeroDim


	} // algo

} // bertini



