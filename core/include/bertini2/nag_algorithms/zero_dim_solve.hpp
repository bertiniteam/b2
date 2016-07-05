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


namespace bertini {

	namespace algorithm {

		template<typename TrackerType, typename EndgameType, typename StartSystemType>
		struct ZeroDim : public Observable<>
		{

			using BaseComplexType = TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    = TrackerTraits<TrackerType>::BaseRealType;
			
			using PrecisionConfig = TrackerTraits<TrackerType>::PrecisionConfig;

			ZeroDim(System const& sys) : target_system_(sys)
			{
				ConsistencyCheck();
			}

			ZeroDim(TrackerType const& t, 
			              EndgameType const& e,
			              System const& s) : target_system_(s), tracker_(t), endgame_(e)
			{
				ConsistencyCheck();
			}



			/**
			\brief Check to ensure that target system is valid for solving.
			*/
			void ConsistencyCheck()
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
				using ComplexFromString = NumTraits<BCT>::FromString;
				using RealFromString = NumTraits<BRT>::FromString;

				target_system_.Homogenize(); // work over projective coordinates
				target_system_.AutoPatch(); // then patch if needed

				start_system_ = StartSystemType(target_system_); // make the start system from the target system.
				start_system_.Homogenize(); //ensure this is homogeneous, too.

				Var t = std::make_shared<Variable>("zero_dim_path_variable"); 

				homotopy_ = (1-t)*target_system_ + std::make_shared<Node>(node::Rational::Rand())*t*start_system_;
				homotopy_.AddPathVariable(t);


				tracker_ = TrackerType(homotopy_);

				tracker_.Setup(DefaultPredictor(),
				              	endgame_.tolerances_.newton_before_endgame, RealFromString("1e5"),
								config::Stepping<BRT>(), config::Newton());

				tracker_.PrecisionSetup(PrecisionConfig(homotopy_));
				
				
				t_start_ = static_cast<BCT>(1)
				t_endgame_boundary_ = ComplexFromString("0.1");
			}

			void Solve()
			{
				auto num_paths_to_track = start_system_.NumStartPoints();
				std::vector<Vec<BCT> > solutions_at_endgame_boundary;
				for (decltype(num_paths_to_track) ii = 0; ii < num_paths_to_track; ++ii)
				{
					DefaultPrecision(ambient_precision);
					homotopy_.precision(ambient_precision);
					auto start_point = start_system_.StartPoint<BCT>(ii);

					Vec<BCT> result;
					SuccessCode tracking_success = tracker.TrackPath(result,t_start,t_endgame_boundary,start_point);

					solutions_at_endgame_boundary.push_back(result);
				}

				MidPathCheck(solutions_at_endgame_boundary);

				tracker.Setup(TestedPredictor,
				              	RealFromString("1e-6"), RealFromString("1e5"),
								stepping_preferences, newton_preferences);


				std::vector<Vec<BCT> > endgame_solutions;
				for (const auto& s : solutions_at_endgame_boundary)
				{
					DefaultPrecision(Precision(s));
					homotopy_.precision(Precision(s));
					SuccessCode endgame_success = endgame_.Run(t_endgame_boundary,s);

					endgame_solutions.push_back(endgame_.FinalApproximation<BCT>())

					if(endgame_success == SuccessCode::Success)
					{

					}
				}

				PostProcessing();
			}


		private:


			void PostProcessing()
			{

			}

		

			PostProcessing<BaseRealType> post_processing_;

			PrecisionConfig precision_config_;
			TrackerType tracker_;
			EndgameType endgame_;

			const System& target_system_;
			StartSystemType start_system_;
			System homotopy_;
		}


	} // algo

} // bertini



