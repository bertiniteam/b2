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

#ifndef BERTINI_TRACKING_CONFIG_HPP
#define BERTINI_TRACKING_CONFIG_HPP


#include <eigen3/Eigen/Dense>
#include "eigen_extensions.hpp"

#include "system.hpp"

namespace bertini
{
	namespace tracking{

		template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;
		template<typename NumType> using Mat = Eigen::Matrix<NumType, Eigen::Dynamic, Eigen::Dynamic>;




		enum class SuccessCode
		{
			Success,
			HigherPrecisionNecessary,
			ReduceStepSize,
			GoingToInfinity,
			FailedToConverge,
			MatrixSolveFailure,
			Failure
		};

		


		namespace config{

			using mpfr_float = boost::multiprecision::mpfr_float;

			enum class Predictor
			{
				Euler
			};

			enum class PrecisionType
			{
				Double,
				FixedMultiple,
				Adaptive
			};







			struct Tolerances
			{
				mpfr_float newton_before_endgame;
				mpfr_float newton_during_endgame;

				mpfr_float final_tolerance;

				mpfr_float path_truncation_threshold;
			};

			struct Stepping
			{
				mpfr_float max_step_size;
				mpfr_float min_step_size;

				mpfr_float step_size_success_factor,
				mpfr_float step_size_fail_factor,

				unsigned consecutive_successful_steps_before_stepsize_increase;
				unsigned consecutive_successful_steps_before_precision_decrease;
				unsigned max_num_steps;
			};

			struct Newton
			{
				unsigned max_num_newton_iterations;
				unsigned min_num_newton_iterations;
			};


			
			struct RenameMe
			{
				unsigned frequency_of_CN_estimation
				
				unsigned sharpendigits;

				mpfr_float function_residual_tolerance;
				mpfr_float ratio_tolerance;
			};


			struct Security
			{
				int level;
				mpfr_float max_norm;
			};






			struct PowerSeries
			{
				mpfr_float SampleFactor;
				unsigned max_cycle_number;
			};

			struct Cauchy
			{
				mpfr_float cutoff_cycle_time;
				mpfr_float cutoff_ratio_time;
			};


			struct TrackBack
			{
				unsigned minimum_cycle;
				bool junk_removal_test;
				unsigned max_depth_LDT;
			};








			struct PostProcessing{
				mpfr_float real_threshold;
				mpfr_float endpoint_finite_threshold;
				mpfr_float final_tol_multiplier;
				mpfr_float final_tol_times_mult;
			};










			struct Regeneration
			{
				bool remove_infinite_endpoints;
				bool higher_dimension_check;

				mpfr_float newton_before_endgame;
				mpfr_float newton_during_endgame;
				mpfr_float final_tolerance;
			};

			





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
				mpfr_float coefficient_bound;  ///< User-defined bound on the sum of the abs vals of the coeffs for any polynomial in the system (for adaptive precision). 
				mpfr_float degree_bound; ///<  User-set bound on degrees of polynomials in the system - tricky to compute for factored polys, subfuncs, etc. (for adaptive precision). 
				mpfr_float epsilon;  ///< Bound on \f$\epsilon\f$ (an error bound).  Used for AMP criteria A, B.
				mpfr_float Phi;  ///< Bound on \f$\Phi\f$ (an error bound).   Used for AMP criteria A, B.
				mpfr_float Psi;  ///< Bound on \f$\Psi\f$ (an error bound).   Used for AMP criterion C.

				int safety_digits_1; ///< User-chosen setting for the number of safety digits used during Criteria A & B.
				int safety_digits_2; ///< User-chosen setting for the number of safety digits used during Criterion C.
				unsigned int maximum_precision; ///< User-chosed setting for the maximum allowable precision.  Paths will die if their precision is requested to be set higher than this threshold.
			
				AdaptiveMultiplePrecisionConfig() : coefficient_bound("1000.0"), degree_bound("5.0"), safety_digits_1(1), safety_digits_2(1), maximum_precision(300) 
				{}

				/**
				 \brief Set epsilon, degree bound, and coefficient bound from system.
				 
				 -Epsilon is set as the square of the number of variables.
				 -Bound on degree is set from a call to System class.  Let this be \f$D\f$  \see System::DegreeBound().
				 -Bound on absolute values of coeffs is set from a call to System class.  Let this be \f$B\f$.  \see System::CoefficientBound().
				*/
				void SetBoundsAndEpsilonFrom(System const& sys)
				{
					epsilon = pow(mpfr_float(sys.NumVariables()),2);
					degree_bound = sys.DegreeBound();
					coefficient_bound = sys.CoefficientBound();
				}
				

				/**
				 Sets values epsilon, Phi, Psi, degree_bound, and coefficient_bound from input system.
				
				 	-Phi becomes \f$D*(D-1)*B\f$.
				 	-Psi is set as \f$D*B\f$.
				*/
				void SetPhiPsiFromBounds()
				{
					Phi = degree_bound*(degree_bound-mpfr_float("1.0"))*coefficient_bound;
				    Psi = degree_bound*coefficient_bound;  //Psi from the AMP paper.
				}

				void SetAMPConfigFrom(System const& sys)
				{
					SetBoundsAndEpsilonFrom(sys);
					SetPhiPsiFromBounds();
				}
			};

			inline
			std::ostream& operator<<(std::ostream & out, AdaptiveMultiplePrecisionConfig const& AMP)
			{
				out << "coefficient_bound: " << AMP.coefficient_bound << "\n";
				out << "degree_bound: " << AMP.degree_bound << "\n";
				out << "epsilon: " << AMP.epsilon << "\n";
				out << "Phi: " << AMP.Phi << "\n";
				out << "Psi: " << AMP.Psi;
				return out;
			}

			
			/**
			\brief Construct a ready-to-go set of AMP settings from a system.

			\see AdaptiveMultiplePrecisionConfig::SetBoundsAndEpsilonFrom
			\see AdaptiveMultiplePrecisionConfig::SetPhiPsiFromBounds
			\see AdaptiveMultiplePrecisionConfig::SetAMPConfigFrom
			*/
			inline
			AdaptiveMultiplePrecisionConfig AMPConfigFrom(System const& sys) 
			{
				AdaptiveMultiplePrecisionConfig AMP;				
				AMP.SetAMPConfigFrom(sys);
			    return AMP;
			}
		}
	}
}


#endif
