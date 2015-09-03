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
#include "eigen_extensions.hpp"


namespace bertini
{
	namespace tracking{

		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;


		enum class SuccessCode
		{
			Success,
			HigherPrecisionNecessary,
			ReduceStepSize,
			GoingToInfinity,
			FailedToConverge,
			MatrixSolveFailure
		};

		


		namespace config{

			enum class PrecisionType
			{
				Double,
				FixedMultiple,
				Adaptive
			};



			namespace general{			
				struct TrackingTolerances
				{
					double before_endgame;
					double during_endgame;
				};
			}

			/**
			Holds the program parameters with respect to Adaptive Multiple Precision.
			
			These criteria are developed in \cite{amp1, amp2}.

			Let:
			\f$J\f$ be the Jacobian matrix of the square system being solved.  
			\f$d\f$ is the latest Newton residual.
			\f$N\f$ is the maximum number of Newton iterations to perform.

			Criterion A:
			\f$ P > \sigma_1 + \log_{10} [ ||J^{-1}|| \epsilon (||J|| + \Phi)   ]  \f$
			
			Criterion B:
			\f$ P > \sigma_1 + D + (\tau + \log_{10} ||d||) / (N-i)  \f$
			where 
			\f$ D = \log_{10} [||J^{-1}||((2 + \epsilon)||J|| + \epsilon \Phi) | 1] \f$

			Criterion C:
			\f$ P > \sigma_2 + \tau + \log_{10}(||J^{-1}|| \Psi + ||z||)  \f$

	
			
			*/
			struct AdaptiveMultiplePrecisionConfig
			{
				double bound_on_abs_vals_of_coeffs;  ///< User-defined bound on the sum of the abs vals of the coeffs for any polynomial in the system (for adaptive precision). 
				double bound_on_degree; ///<  User-set bound on degrees of polynomials in the system - tricky to compute for factored polys, subfuncs, etc. (for adaptive precision). 
				double epsilon;  ///< Bound on \f$\epsilon\f$ (an error bound).  Used for AMP criteria A, B.
				double Phi;  ///< Bound on \f$\Phi\f$ (an error bound).   Used for AMP criteria A, B.
				double Psi;  ///< Bound on \f$\Psi\f$ (an error bound).   Used for AMP criterion C.

				int safety_digits_1;
				int safety_digits_2;
				int maximum_precision;
			} ;
		}
	}
}


#endif
