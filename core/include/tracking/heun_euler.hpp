//This file is part of Bertini 2.0.
//
//heun_euler.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//heun_euler.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with heun_euler.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  heun_euler.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Fall 2015

#ifndef BERTINI_HEUN_EULER_PREDICT_HPP
#define BERTINI_HEUN_EULER_PREDICT_HPP


#include "tracking/amp_criteria.hpp"
#include "tracking/tracking_config.hpp"

#include "system.hpp"
#include <eigen3/Eigen/LU>

#include <boost/type_index.hpp>


namespace bertini{
	namespace tracking{
		namespace predict{



			// struct HeunEuler{

			// 	template <typename ComplexType, typename RealType>
			// 	SuccessCode operator(Vec<ComplexType> & next_space,
			// 		                  RealType & error_estimate,
			// 		                  RealType & size_proportion,
			// 		                  RealType & norm_J,
			// 		                  RealType & norm_J_inverse,
			// 		               System & S,
			// 		               Vec<ComplexType> const& current_space, ComplexType current_time, 
			// 		               ComplexType const& delta_t,
			// 		               RealType & condition_number_estimate,
			// 		               unsigned & num_steps_since_last_condition_number_computation, 
			// 		               unsigned frequency_of_CN_estimation, 
			// 		               RealType const& tracking_tolerance,
			// 		               config::AdaptiveMultiplePrecisionConfig const& AMP_config)

			template <typename ComplexType, typename RealType>
			SuccessCode HeunEuler(Vec<ComplexType> & next_space,
				                  RealType & error_estimate,
				                  RealType & size_proportion,
				                  RealType & norm_J,
				                  RealType & norm_J_inverse,
				               System & S,
				               Vec<ComplexType> const& current_space, ComplexType current_time, 
				               ComplexType const& delta_t,
				               RealType & condition_number_estimate,
				               unsigned & num_steps_since_last_condition_number_computation, 
				               unsigned frequency_of_CN_estimation, 
				               RealType const& tracking_tolerance,
				               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				Mat<ComplexType> dh_dx = S.Jacobian(current_space, current_time);
				auto LU = dh_dx.lu(); // we keep the LU here because need to estimate the condition number of J^-1
				
				if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
					return SuccessCode::MatrixSolveFailure;

				auto delta_x_1 = LU.solve(-S.TimeDerivative(current_space, current_time)); 

				norm_J = dh_dx.norm();
				norm_J_inverse = LU.solve(Vec<ComplexType>::Random(S.NumVariables())).norm();

				condition_number_estimate = norm_J*norm_J_inverse;






				Vec<ComplexType> second_sample_point = current_space + delta_x_1*delta_t;
				ComplexType second_time = current_time + delta_t;

				LU = S.Jacobian(second_sample_point, second_time).lu(); // we keep the LU here because need to estimate the condition number of J^-1
				
				if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
					return SuccessCode::MatrixSolveFailure;

				auto delta_x_2 = LU.solve(-S.TimeDerivative(second_sample_point, second_time)); 

				auto delta_x = (delta_x_1+delta_x_2)/2;

				auto Y = (delta_x_1 - delta_x_2)/2;

				error_estimate = Y.norm();

				if (error_estimate==RealType(0) 
				    || 
				    -log10(error_estimate) >= NumTraits<RealType>::NumFuzzyDigits()
				    ||
				    -log10(abs2(delta_t)) >= 2*NumTraits<RealType>::NumFuzzyDigits()
				    )
				{
					size_proportion = RealType(1);
				}
				else
				{
					size_proportion = pow(10, log10(error_estimate) - log10(abs2(delta_t))); 
				}

				return SuccessCode::Success;
			}

		} //re: namespace predict
	} // re: namespace tracking
} // re: namespace bertini
#endif


