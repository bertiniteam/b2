//This file is part of Bertini 2.0.
//
//correct.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//correct.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with correct.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  correct.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015


/**
\file correct.hpp 

\brief Contains wrappers for the Newton correct functions.
*/

#ifndef BERTINI_CORRECT_HPP
#define BERTINI_CORRECT_HPP

#include "tracking/newton_correct.hpp"

namespace bertini{
	namespace tracking{
			
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
		SuccessCode Correct(Vec<ComplexType> & next_space,
					               System const& S,
					               Vec<ComplexType> const& current_space, 
					               ComplexType const& current_time, 
					               RealType tracking_tolerance,
					               RealType path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations)
		{
			return correct::NewtonLoop(next_space,
						               S,
						               current_space, 
						               current_time, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               min_num_newton_iterations,
						               max_num_newton_iterations);
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
		SuccessCode Correct(Vec<ComplexType> & next_space,
					               System const& S,
					               Vec<ComplexType> const& current_space, 
					               ComplexType const& current_time, 
					               RealType tracking_tolerance,
					               RealType path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			return correct::NewtonLoop(next_space,
						               S,
						               current_space, 
						               current_time, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               min_num_newton_iterations,
						               max_num_newton_iterations,
						               AMP_config);
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
		SuccessCode Correct(Vec<ComplexType> & next_space,
		                    RealType & norm_delta_z,
			                       RealType & norm_J,
			                       RealType & norm_J_inverse,
			                       RealType & condition_number_estimate,
					               System const& S,
					               Vec<ComplexType> const& current_space, 
					               ComplexType const& current_time, 
					               RealType tracking_tolerance,
					               RealType path_truncation_threshold,
					               unsigned min_num_newton_iterations,
					               unsigned max_num_newton_iterations,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
		{
			return correct::NewtonLoop(next_space,
			                           norm_delta_z,
			                           norm_J,
			                           norm_J_inverse,
			                           condition_number_estimate,
						               S,
						               current_space, 
						               current_time, 
						               tracking_tolerance,
						               path_truncation_threshold,
						               min_num_newton_iterations,
						               max_num_newton_iterations,
						               AMP_config);
		}


	}
	
}


#endif


