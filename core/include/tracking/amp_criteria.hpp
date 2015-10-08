//This file is part of Bertini 2.0.
//
//amp_criteria.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_criteria.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_criteria.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  amp_criteria.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef BERTINI_AMP_CRITERIA_HPP
#define BERTINI_AMP_CRITERIA_HPP

#include <boost/multiprecision/mpfr.hpp>

#include "tracking/tracking_config.hpp"

namespace bertini{
	namespace tracking{
		namespace amp{



			using AdaptiveMultiplePrecisionConfig = config::AdaptiveMultiplePrecisionConfig;


			/**
			Check AMP Criterion A, from \cite{amp1, amp2}.

			True means the check passed, and the precision is all good.  False means something's gotta change, stepsize or precision.
			
			\param norm_J The matrix norm of the Jacoabian matrix
			\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian matrix.
			\param AMP_config The settings for adaptive multiple precision.

			\return True if criteria satisfied, false if violated and precision or step length should be adjusted.
			*/
			template<typename RealType>
			bool CriterionA(RealType const& norm_J, RealType const& norm_J_inverse, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				// std::cout << "criterion A, lhs: " << NumTraits<RealType>::NumDigits() << "\n";
				// std::cout << "criterion A, rhs: " << AMP_config.safety_digits_1 + log10(norm_J_inverse * AMP_config.epsilon * (norm_J + AMP_config.Phi)) << "\n";
				return NumTraits<RealType>::NumDigits()  > AMP_config.safety_digits_1 + log10(norm_J_inverse * RealType(AMP_config.epsilon) * (norm_J + RealType(AMP_config.Phi) ) );
			}
			


			/**
			\brief Compute the expression \f$D\f$ from the AMP papers \cite{AMP1,AMP2}.
			*/
			template <typename RealType>
			RealType D(RealType const& norm_J, RealType const& norm_J_inverse, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				return log10(norm_J_inverse*( (RealType(2)+RealType(AMP_config.epsilon))*norm_J + RealType(AMP_config.epsilon)*RealType(AMP_config.Phi)) + RealType(1));
			}

			/**
			Evaluate the right hand side of the inequality of Criterion B, from \cite{amp1, amp2}.
			
			\param norm_J The matrix norm of the Jacoabian matrix
			\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian matrix.
			\param num_newton_iterations_remaining The number of iterations which have yet to perform.
			\param tracking_tolerance The tightness to which the path should be tracked.  This is the raw tracking tolerance (for now)
			\param latest_newton_residual The norm of the length of the most recent Newton step.
			\param AMP_config The settings for adaptive multiple precision.

			\return The value of the right hand side of Criterion B
			*/
			template<typename RealType>
			RealType CriterionBRHS(RealType const& norm_J, RealType const& norm_J_inverse, unsigned num_newton_iterations_remaining, RealType const& tracking_tolerance, RealType const& norm_of_latest_newton_residual, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				// std::cout << norm_J << "\n";
				// std::cout << norm_J_inverse << "\n";
				// std::cout << num_newton_iterations_remaining << "\n";
				// std::cout << tracking_tolerance << "\n";
				// std::cout << norm_of_latest_newton_residual << "\n";
				// std::cout << AMP_config << "\n";
				
				// std::cout << D << "\n";
				// std::cout << "criterion B, lhs: " << NumTraits<RealType>::NumDigits() << "\n";
				// std::cout << "criterion B, rhs: " << RealType(AMP_config.safety_digits_1 + D + (-log10(tracking_tolerance) + log10(norm_of_latest_newton_residual)) / (num_newton_iterations_remaining)) << "\n";
				return AMP_config.safety_digits_1 + D(norm_J, norm_J_inverse, AMP_config) + (-log10(tracking_tolerance) + log10(norm_of_latest_newton_residual)) / (num_newton_iterations_remaining);
			}


			/**
			Check AMP Criterion B, from \cite{amp1, amp2}.

			True means the check passed, and the precision is all good.  False means something's gotta change, stepsize or precision.
			
			\param norm_J The matrix norm of the Jacoabian matrix
			\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian matrix.
			\param num_newton_iterations_remaining The number of iterations which have yet to perform.
			\param tracking_tolerance The tightness to which the path should be tracked.  This is the raw tracking tolerance (for now)
			\param latest_newton_residual The norm of the length of the most recent Newton step.
			\param AMP_config The settings for adaptive multiple precision.

			\return True if criteria satisfied, false if violated and precision or step length should be adjusted.
			*/
			template<typename RealType>
			bool CriterionB(RealType const& norm_J, RealType const& norm_J_inverse, unsigned num_newton_iterations_remaining, RealType const& tracking_tolerance, RealType const& norm_of_latest_newton_residual, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				return NumTraits<RealType>::NumDigits() > CriterionBRHS(norm_J, norm_J_inverse, num_newton_iterations_remaining, tracking_tolerance,  norm_of_latest_newton_residual, AMP_config);
			}


			



			/**
			Evaluate the right hand side of Criterion C, from \cite{amp1, amp2}.

			\param norm_J The matrix norm of the Jacoabian matrix
			\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian matrix.
			\param tracking_tolerance The tightness to which the path should be tracked.  This is the raw tracking tolerance
			\param z The current space point?  TODO: check this. 
			\param AMP_config The settings for adaptive multiple precision.

			\return The value of the right hand side of the inequality from Criterion C.
			*/
			template<typename ComplexType ,typename RealType>
			RealType CriterionCRHS(RealType const& norm_J_inverse, Vec<ComplexType> const& z, RealType tracking_tolerance, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");
				// std::cout << "criterion C, lhs: " << NumTraits<RealType>::NumDigits() << "\n";
				// std::cout << "criterion C, rhs: " << AMP_config.safety_digits_2 + -log10(tracking_tolerance) + log10(norm_J_inverse*AMP_config.Psi + z.norm()) << "\n";
				return AMP_config.safety_digits_2 + -log10(tracking_tolerance) + log10(norm_J_inverse*RealType(AMP_config.Psi) + z.norm());
			}


			/**
			Check AMP Criterion C, from \cite{amp1, amp2}.

			True means the check passed, and the precision is all good.  False means something's gotta change, stepsize or precision.
			
			\param norm_J The matrix norm of the Jacoabian matrix
			\param norm_J_inverse An estimate on the norm of the inverse of the Jacobian matrix.
			\param tracking_tolerance The tightness to which the path should be tracked.  This is the raw tracking tolerance
			\param z The current space point?  TODO: check this. 
			\param AMP_config The settings for adaptive multiple precision.
			*/
			template<typename ComplexType ,typename RealType>
			bool CriterionC(RealType const& norm_J_inverse, Vec<ComplexType> const& z, RealType tracking_tolerance, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				static_assert(std::is_same<typename Eigen::NumTraits<RealType>::Real, typename Eigen::NumTraits<ComplexType>::Real>::value,"underlying complex type and the type for comparisons must match");

				return NumTraits<RealType>::NumDigits() > CriterionCRHS(norm_J_inverse, z, tracking_tolerance, AMP_config);
			}
		}
	}
	
}


#endif

