//This file is part of Bertini 2.0.
//
//euler.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//euler.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with euler.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  euler.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef euler_hpp
#define euler_hpp

#include "tracking/amp_criteria.hpp"
#include "tracking/tracking_config.hpp"

#include "system.hpp"
#include <eigen3/Eigen/LU>

#include <boost/type_index.hpp>

namespace bertini{
	namespace tracking{





		namespace predict{

			using PrecisionType = config::PrecisionType;
			/**
			Perform an euler-prediction step

			\param predictor_choice The enum class selecting the predictor to be used.
			\param next_space The computed prediction.
			\param sys The system being solved.
			\param current_space The current space variable vector.
			\param current_time The current time.
			\param delta_t The size of the time step.
			\param condition_number_estimate The computed estimate of the condition number of the Jacobian.
			\param num_steps_since_last_condition_number_computation.  Updated in this function.
			\param frequency_of_CN_estimation How many steps to take between condition number estimates.
			\param prec_type The operating precision type.  
			\param tracking_tolerance How tightly to track the path.
			\param AMP_config The settings for adaptive multiple precision.

			\tparam ComplexType The number type for evaluation.
			\tparam RealType The number type for comparisons.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode Euler(Vec<ComplexType> & next_space,
					               System & S,
					               Vec<ComplexType> const& current_space, ComplexType current_time, 
					               ComplexType const& delta_t,
					               RealType & condition_number_estimate,
					               unsigned & num_steps_since_last_condition_number_computation, 
					               unsigned frequency_of_CN_estimation, PrecisionType PrecType, 
					               RealType const& tracking_tolerance,
					               config::AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

				// std::cout << "current_time = " << current_time << "\n";
				// std::cout << "current_space = " << current_space << "\n";

				Vec<ComplexType> dh_dt = -S.TimeDerivative(current_space, current_time);
				Mat<ComplexType> dh_dx = S.Jacobian(current_space, current_time); // this will complain (throw) if the system does not depend on time.

				// std::cout << "size of dh_dt = \n" << dh_dt.size() << "\n";
				// std::cout << "type of dh_dt = \n" << boost::typeindex::type_id_with_cvr<decltype(dh_dt)>().pretty_name() << "\n";
				// std::cout << "dh_dt = \n" << dh_dt << "\n";
				// std::cout << "dh_dx = \n" << dh_dx << "\n";

				// solve delta_x = (dH/dx)^(-1)*Y
				// Y = dH/dt

				auto LU = dh_dx.lu(); // we keep the LU here because may need to estimate the condition number of J^-1
				
				if (LUPartialPivotDecompositionSuccessful(LU.matrixLU())!=MatrixSuccessCode::Success)
					{
						if (PrecType==PrecisionType::Adaptive)
							return SuccessCode::HigherPrecisionNecessary;
						else
							return SuccessCode::MatrixSolveFailure;
					}

				auto delta_x = LU.solve(dh_dt); 

				// std::cout << "euler delta_x = \n" << delta_x << std::endl;

				if (PrecType==PrecisionType::Adaptive)
				{
					Vec<ComplexType> randy = Vec<ComplexType>::Random(S.NumVariables());
					Vec<ComplexType> temp_soln = LU.solve(randy);
					
					RealType norm_J = dh_dx.norm();
					RealType norm_J_inverse = temp_soln.norm();


					if (num_steps_since_last_condition_number_computation >= frequency_of_CN_estimation)
					{
						// TODO: this esimate may be wrong
						condition_number_estimate = norm_J * norm_J_inverse;
						num_steps_since_last_condition_number_computation = 0; // reset the counter to 0
					}
					else // no need to compute the condition number
						num_steps_since_last_condition_number_computation++;


					// std::cout << "norm_J: " << norm_J << " norm_J_inverse: " << norm_J_inverse << "\n";
					// std::cout << "condition_number_estimate: " << condition_number_estimate << "\n";

					if (!amp::CriterionA(norm_J, norm_J_inverse, AMP_config)) // AMP_criterion_A != ok
						return SuccessCode::HigherPrecisionNecessary;
					else if (!amp::CriterionC(norm_J_inverse, current_space, tracking_tolerance, AMP_config)) // AMP_criterion_C != ok
						return SuccessCode::HigherPrecisionNecessary;

				}

				// std::cout << "euler delta_t * delta_x = \n" << delta_x * delta_t << std::endl;

				next_space = current_space + delta_x * delta_t;

				return SuccessCode::Success;
			}
		}
	}
}



#endif

