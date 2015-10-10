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

#ifndef BERTINI_NEWTON_CORRECT_HPP
#define BERTINI_NEWTON_CORRECT_HPP


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
			template <typename ComplexType, typename RealType>
			SuccessCode NewtonLoop(Vec<ComplexType> & next_space,
					               System & S,
					               Vec<ComplexType> const& current_space, // pass by value to get a copy of it
					               ComplexType const& current_time, 
					               RealType const& tracking_tolerance,
					               RealType const& path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(max_num_newton_iterations >= min_num_newton_iterations && "max number newton iterations must be at least the min.");
				#endif


				next_space = current_space;
				for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
				{
					// std::cout << ii << "th iteration\n";
					//TODO: wrap these into a single line.
					auto f = S.Eval(next_space, current_time);
					auto J = S.Jacobian(next_space, current_time);
					auto LU = J.lu();


					if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;


					auto delta_z = LU.solve(-f);
					// std::cout << "correct delta_z = \n" << delta_z << std::endl;
					next_space += delta_z;

					if ( delta_z.norm() < tracking_tolerance && ii >= (min_num_newton_iterations-1) )
						return SuccessCode::Success;


					auto norm_J_inverse = LU.solve(Vec<ComplexType>::Random(S.NumVariables())).norm();
					if (!amp::CriterionB(J.norm(), norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, delta_z.norm(), AMP_config))
						return SuccessCode::HigherPrecisionNecessary;

					if (!amp::CriterionC(norm_J_inverse, next_space, tracking_tolerance, AMP_config))
						return SuccessCode::HigherPrecisionNecessary;

					if (S.DehomogenizePoint(next_space).norm() > path_truncation_threshold)
						return SuccessCode::GoingToInfinity;

				}

				return SuccessCode::FailedToConverge;
			}
			
		}
		
	}
	
}

#endif
