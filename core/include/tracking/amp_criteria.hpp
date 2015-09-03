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

#ifndef amp_criteria_hpp
#define amp_criteria_hpp

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
			template<typename T>
			bool CriterionA(T const& norm_J, T const& norm_J_inverse, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				return boost::multiprecision::mpfr_float::default_precision() > AMP_config.safety_digits_1 + log10(norm_J_inverse * AMP_config.epsilon * (norm_J + AMP_config.Phi));
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
			template<typename T>
			bool CriterionB(T const& norm_J, T const& norm_J_inverse, unsigned num_newton_iterations_remaining, double tracking_tolerance, Vec<T> const& latest_newton_residual, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				auto D = log10(norm_J_inverse*( (2+AMP_config.epsilon)*norm_J + AMP_config.epsilon*AMP_config.Phi) + 1);
				return boost::multiprecision::mpfr_float::default_precision() > AMP_config.safety_digits_1 + D + (-log10(tracking_tolerance) + log10(norm(latest_newton_residual))) / (num_newton_iterations_remaining);
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
			template<typename T>
			bool CriterionC(T const& norm_J_inverse, Vec<T> const& z, double tracking_tolerance, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				return boost::multiprecision::mpfr_float::default_precision() > AMP_config.safety_digits_2 + -log10(tracking_tolerance) + log10(norm_J_inverse*AMP_config.Psi + norm(z));
			}
		}
	}
	
}


#endif

