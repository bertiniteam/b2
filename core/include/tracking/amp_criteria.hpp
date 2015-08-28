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

#include "tracking/tracking_config.hpp"

namespace bertini{
	namespace tracking{
		namespace amp{

			inline
			template<typename T>
			bool CriterionA(T norm_J, T norm_J_inverse, AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				return boost::multiprecision::mpfr_float::default_precision() > AMP_config.safety_digits_1 + log10(norm_J_inverse * AMP_config.epsilon * (norm_J + AMP_config.Phi));
			}
			
		}
	}
	
}


#endif

