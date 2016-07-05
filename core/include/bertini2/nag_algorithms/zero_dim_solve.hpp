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

namespace bertini {

	namespace algorithm {

		template<typename TrackerType, typename EndgameType, typename StartSystemType>
		struct ZeroDim : public Observable<>
		{
			BERTINI_DEFAULT_VISITABLE();

			using BaseComplexType = typename tracking::TrackerTraits<TrackerType>::BaseComplexType;
			using BaseRealType    = typename tracking::TrackerTraits<TrackerType>::BaseRealType;
			
			using PrecisionConfig = typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;

			ZeroDim(System const& sys) : target_system_(sys), tracker_(sys), endgame_(tracker_)
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
				using Variable = node::Variable;
				using Var = std::shared_ptr<Variable>;

				target_system_.Homogenize(); // work over projective coordinates
				target_system_.AutoPatch(); // then patch if needed

				start_system_ = StartSystemType(target_system_); // make the start system from the target system.
				start_system_.Homogenize(); //ensure this is homogeneous, too.

				Var t = std::make_shared<Variable>("zero_dim_path_variable"); 

				homotopy_ = (1-t)*target_system_ + std::make_shared<node::Rational>(node::Rational::Rand())*t*start_system_;
				homotopy_.AddPathVariable(t);

				auto tolerances = algorithm::config::Tolerances<BaseRealType>();

				tracker_ = TrackerType(homotopy_);

				tracker_.Setup(tracking::predict::DefaultPredictor(),
				              	tolerances.newton_before_endgame, NumTraits<BaseRealType>::FromString("1e5"),
								tracking::config::Stepping<BaseRealType>(), tracking::config::Newton());

				tracker_.PrecisionSetup(PrecisionConfig(homotopy_));
				
				
				t_start_ = static_cast<BaseComplexType>(1);
				t_endgame_boundary_ = NumTraits<BaseComplexType>::FromString("0.1");
				t_end_ = static_cast<BaseComplexType>(0);

				ambient_precision_ = DoublePrecision();
			}

			void Solve()
			{
				using SuccessCode = tracking::SuccessCode;


				auto num_paths_to_track = start_system_.NumStartPoints();
				std::vector<Vec<BaseComplexType> > solutions_at_endgame_boundary;
				for (decltype(num_paths_to_track) ii = 0; ii < num_paths_to_track; ++ii)
				{
					DefaultPrecision(ambient_precision_);
					homotopy_.precision(ambient_precision_);
					auto start_point = start_system_.template StartPoint<BaseComplexType>(ii);

					Vec<BaseComplexType> result;
					SuccessCode tracking_success = tracker_.TrackPath(result,t_start_,t_endgame_boundary_,start_point);

					solutions_at_endgame_boundary.push_back(result);
				}

				MidPathCheck(solutions_at_endgame_boundary);

				tracker_.SetTrackingTolerance(NumTraits<BaseRealType>::FromString("1e-6"));


				std::vector<Vec<BaseComplexType> > endgame_solutions;
				for (const auto& s : solutions_at_endgame_boundary)
				{
					DefaultPrecision(Precision(s));
					homotopy_.precision(Precision(s));
					SuccessCode endgame_success = endgame_.Run(t_endgame_boundary_,s);

					endgame_solutions.push_back(endgame_.template FinalApproximation<BaseComplexType>());

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

			template <typename T>
			void MidPathCheck(T const& s)
			{

			}


			
			config::PostProcessing<BaseRealType> post_processing_;

			PrecisionConfig precision_config_;
			TrackerType tracker_;
			EndgameType endgame_;

			System target_system_;
			StartSystemType start_system_;
			System homotopy_;

			BaseComplexType t_start_, t_endgame_boundary_, t_end_;

			unsigned ambient_precision_ = DoublePrecision();
		}; // struct ZeroDim


	} // algo

} // bertini



