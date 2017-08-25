//This file is part of Bertini 2.
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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


/**
\file newton_correct.hpp 

\brief Implements Newton's method for square systems
*/


#ifndef BERTINI_NEWTON_CORRECT_HPP
#define BERTINI_NEWTON_CORRECT_HPP


#include "bertini2/trackers/amp_criteria.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/system.hpp"

namespace bertini{
	namespace tracking{
		namespace correct{

			

			/**
			\brief Run Newton's method in fixed precision.

			Run Newton's method until it converges (\f$\Delta z\f$ < tol), or the next point's norm exceeds the path truncation threshold.
	
			\return The SuccessCode indicating what happened.

			\tparam ComplexType The complex type for arithmetic
			\tparam RealType The underlying real number type, used for comparitors.

			\param[out] next_space The computed next space point.
			\param S The system we are tracking on.
			\param current_space The base point for newton correcting.
			\param current_time The current time value.  Note it is complex.
			\param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
			\param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
			\param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
			\param max_num_newton_iterations The maximum number of iterations to run Newton's method for.

			*/
			template <typename ComplexType, typename RealType>
			SuccessCode NewtonLoop(Vec<ComplexType> & next_space,
					               System const& S,
					               Vec<ComplexType> const& current_space, // pass by value to get a copy of it
					               ComplexType const& current_time, 
					               RealType const& tracking_tolerance,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations)
			{
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(max_num_newton_iterations >= min_num_newton_iterations && "max number newton iterations must be at least the min.");
				#endif


				next_space = current_space;
				for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
				{
					auto f = S.Eval(next_space, current_time);
					auto J = S.Jacobian(next_space, current_time);
					auto LU = J.lu();

					if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;

					auto delta_z = LU.solve(-f);
					next_space += delta_z;

					if ( (delta_z.norm() < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
						return SuccessCode::Success;
				}

				return SuccessCode::FailedToConverge;
			}
			



			/**
			\brief Run Newton's method in multiple precision.

			Run Newton's method until it converges (\f$\Delta z\f$ < tol), an AMP criterion (B or C) is violated, or the next point's norm exceeds the path truncation threshold.

			\tparam ComplexType The complex type for arithmetic
			\tparam RealType The underlying real number type, used for comparitors.

			\param[out] next_space The computed next space point.
			\param S The system we are tracking on.
			\param current_space The base point for newton correcting.
			\param current_time The current time value.  Note it is complex.
			\param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
			\param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
			\param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
			\param max_num_newton_iterations The maximum number of iterations to run Newton's method for.
			\param AMP_config Adaptive multiple precision settings.  Using this argument is how Bertini2 knows you want to use adaptive precision.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode NewtonLoop(Vec<ComplexType> & next_space,
					               System const& S,
					               Vec<ComplexType> const& current_space, // pass by value to get a copy of it
					               ComplexType const& current_time, 
					               RealType const& tracking_tolerance,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(max_num_newton_iterations >= min_num_newton_iterations && "max number newton iterations must be at least the min.");
				#endif


				next_space = current_space;
				for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
				{
					//TODO: wrap these into a single line.
					auto f = S.Eval(next_space, current_time);
					auto J = S.Jacobian(next_space, current_time);
					auto LU = J.lu();

					if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;


					auto delta_z = LU.solve(-f);
					next_space += delta_z;

					if ( (delta_z.norm() < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
						return SuccessCode::Success;

					auto norm_J_inverse = LU.solve(RandomOfUnits<ComplexType>(S.NumVariables())).norm();
					if (!amp::CriterionB(J.norm(), norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, delta_z.norm(), AMP_config))
						return SuccessCode::HigherPrecisionNecessary;

					if (!amp::CriterionC(norm_J_inverse, next_space, tracking_tolerance, AMP_config))
						return SuccessCode::HigherPrecisionNecessary;
				}

				return SuccessCode::FailedToConverge;
			}
			






			/**
			\brief Run Newton's method in multiple precision.

			Run Newton's method until it converges (\f$\Delta z\f$ < tol), an AMP criterion (B or C) is violated, or the next point's norm exceeds the path truncation threshold.

			\tparam ComplexType The complex type for arithmetic
			\tparam RealType The underlying real number type, used for comparitors.

			\param[out] next_space The computed next space point.
			\param[out] norm_delta_z The norm of the last step size.
			\param[out] norm_J The matrix norm of the Jacobian matrix.
			\param[out] norm_J_inverse The matrix norm of the inverse of the Jacobian matrix.
			\param[out] condition_number_estimate An estimate on the condition number.
			\param S The system we are tracking on.
			\param current_space The base point for newton correcting.
			\param current_time The current time value.  Note it is complex.
			\param tracking_tolerance The upper threshold for step size.  Must iterate correcting until the corrector step is less than this threshold in length.
			\param path_truncation_threshold Correcting stops the the norm of the current solution exceeds this number.
			\param min_num_newton_iterations The corrector must take at least this many steps.  This should be at least 1.
			\param max_num_newton_iterations The maximum number of iterations to run Newton's method for.
			\param AMP_config Adaptive multiple precision settings.  Using this argument is how Bertini2 knows you want to use adaptive precision.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode NewtonLoop(Vec<ComplexType> & next_space,
			                       RealType & norm_delta_z,
			                       RealType & norm_J,
			                       RealType & norm_J_inverse,
			                       RealType & condition_number_estimate,
					               System const& S,
					               Vec<ComplexType> const& current_space, // pass by value to get a copy of it
					               ComplexType const& current_time, 
					               RealType const& tracking_tolerance,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				#ifndef BERTINI_DISABLE_ASSERTS
				assert(max_num_newton_iterations >= min_num_newton_iterations && "max number newton iterations must be at least the min.");
				#endif


				next_space = current_space;
				for (unsigned ii = 0; ii < max_num_newton_iterations; ++ii)
				{
					//TODO: wrap these into a single line.
					auto f = S.Eval(next_space, current_time);
					auto J = S.Jacobian(next_space, current_time);
					auto LU = J.lu();


					if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
						return SuccessCode::MatrixSolveFailure;


					auto delta_z = LU.solve(-f);
					next_space += delta_z;


					norm_delta_z = delta_z.norm();
					norm_J = J.norm();
					norm_J_inverse = LU.solve(RandomOfUnits<ComplexType>(S.NumVariables())).norm();
					condition_number_estimate = norm_J*norm_J_inverse;



					if ( (norm_delta_z < tracking_tolerance) && (ii >= (min_num_newton_iterations-1)) )
						return SuccessCode::Success;

					if (!amp::CriterionB(norm_J, norm_J_inverse, max_num_newton_iterations - ii, tracking_tolerance, norm_delta_z, AMP_config))
						return SuccessCode::HigherPrecisionNecessary;

					if (!amp::CriterionC(norm_J_inverse, next_space, tracking_tolerance, AMP_config))
						return SuccessCode::HigherPrecisionNecessary;
				}

				return SuccessCode::FailedToConverge;
			}
		}
		
	}
	
}

#endif
