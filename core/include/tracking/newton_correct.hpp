//This file is part of Bertini 2.0.
//
//newton_correct.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//newton_correct.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with newton_correct.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  newton_correct.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef newton_correct_hpp
#define newton_correct_hpp


#include "tracking/amp_criteria.hpp"
#include "tracking/tracking_config.hpp"
#include "system.hpp"

namespace bertini{
	namespace tracking{
		namespace correct{

			using PrecisionType = config::PrecisionType;
			
			/**
			\brief Run Newton's method.

			Run Newton's method until it converges (\f$\Delta z\f$ < tol), an AMP criterion (B or C) is violated, or the next point's norm exceeds the path truncation threshold.
			*/
			template <typename NumType>
			SuccessCode NewtonLoop(Vec<NumType> & next_space,
					               System & S,
					               Vec<NumType> const& current_space, // pass by value to get a copy of it
					               NumType const& current_time, 
					               PrecisionType PrecType, 
					               NumType tracking_tolerance,
					               NumType path_truncation_threshold,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				next_space = current_space;
				for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
				{
					//TODO: wrap these into a single line.
					auto f = S.Eval(current_space, current_time);
					auto J = S.Jacobian(current_space, current_time);
					auto LU = J.lu();


					if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
					{
						if (PrecType==PrecisionType::Adaptive)
							return SuccessCode::HigherPrecisionNecessary;
						else
							return SuccessCode::MatrixSolveFailure;
					}


					auto delta_z = LU.solve(-f);

					if (norm(delta_z) < tracking_tolerance)
						return SuccessCode::Success;

					if (PrecType==PrecisionType::Adaptive)
					{
						auto norm_J_inverse = norm(LU.solve(Vec<NumType>::Random(S.NumVariables())));
						if (!CriterionB(norm(J), norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, delta_z, AMP_config))
							return SuccessCode::HigherPrecisionNecessary;

						if (!CriterionC(norm_J_inverse, current_space, tracking_tolerance, AMP_config))
							return SuccessCode::HigherPrecisionNecessary;

					}

					if (norm(S.DehomogenizePoint(current_space)) > path_truncation_threshold)
					{
						return SuccessCode::GoingToInfinity;
					}

					next_space+=delta_z;
				}

				return SuccessCode::FailedToConverge;
			}
			
		}
		
	}
	
}

#endif
