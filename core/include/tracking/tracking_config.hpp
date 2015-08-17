//This file is part of Bertini 2.0.
//
//tracking_config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking_config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking_config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  tracking_config.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Summer 2015

#ifndef tracking_config_hpp
#define tracking_config_hpp

#include <eigen3/Eigen/Dense>

namespace bertini
{
	namespace tracking{

		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;


		enum class SuccessCode
		{
			Success,
			HigherPrecisionNecessary,
			Failure
		};

		enum class PrecisionType
		{
			Double,
			FixedMultiple,
			Adaptive
		};
	}
}


#endif
