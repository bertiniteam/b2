//This file is part of Bertini 2.
//
//tracking/tracking_config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//tracking/tracking_config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/tracking_config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University

#ifndef BERTINI_TRACKING_CONFIG_HPP
#define BERTINI_TRACKING_CONFIG_HPP

/**
\file tracking_config.hpp

\brief Configs and settings for tracking.
*/
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"

#include "bertini2/system.hpp"

namespace bertini
{
	namespace tracking{

		// aliases for the types used to contain space and time samples, for the endgames.
		template<typename T> using SampCont = std::deque<Vec<T> >;
		template<typename T> using TimeCont = std::deque<T>;
		
		
		enum class PrecisionType
		{
			Fixed,
			Adaptive
		};

		// forward declarations
		template<class D>
		class FixedPrecisionTracker;
		class MultiplePrecisionTracker;
		class DoublePrecisionTracker;

		namespace config{
			template<typename T>
			struct FixedPrecisionConfig;
			struct AdaptiveMultiplePrecisionConfig;
		}
		// end forward declarations




		// now for the TrackerTraits structs, which enable lookup of correct settings objects and types, etc.
		template<class T>
		struct TrackerTraits
		{};



		template<>
		struct TrackerTraits<DoublePrecisionTracker>
		{
			using BaseComplexType = dbl;
			using BaseRealType = double;
			using EventEmitterType = FixedPrecisionTracker<DoublePrecisionTracker>;
			using PrecisionConfig = config::FixedPrecisionConfig<BaseRealType>;
			enum {
				IsFixedPrec = 1,
				IsAdaptivePrec = 0
			};
		};


		template<>
		struct TrackerTraits<MultiplePrecisionTracker>
		{
			using BaseComplexType = mpfr;
			using BaseRealType = mpfr_float;
			using EventEmitterType = FixedPrecisionTracker<MultiplePrecisionTracker>;
			using PrecisionConfig = config::FixedPrecisionConfig<BaseRealType>;

			enum {
				IsFixedPrec = 1,
				IsAdaptivePrec = 0
			};
		};

		class AMPTracker;
		template<>
		struct TrackerTraits<AMPTracker>
		{
			using BaseComplexType = mpfr;
			using BaseRealType = mpfr_float;
			using EventEmitterType = AMPTracker;
			using PrecisionConfig = config::AdaptiveMultiplePrecisionConfig;

			enum {
				IsFixedPrec = 0,
				IsAdaptivePrec = 1
			};
		};


		

		template<class D>
		struct TrackerTraits<FixedPrecisionTracker<D> >
		{
			using BaseComplexType = typename TrackerTraits<D>::BaseComplexType;
			using BaseRealType = typename TrackerTraits<D>::BaseRealType;
			using EventEmitterType = typename TrackerTraits<D>::EventEmitterType;
			using PrecisionConfig = typename TrackerTraits<D>::PrecisionConfig;

			enum {
				IsFixedPrec = 1,
				IsAdaptivePrec = 0
			};
		};


		namespace endgame{
		// some forward declarations
		template<typename Tracker, typename Enable = void>
		class FixedPrecPowerSeriesEndgame;

		template<typename Tracker, typename Enable = void>
		class FixedPrecCauchyEndgame;

		class AMPPowerSeriesEndgame;
		class AMPCauchyEndgame;
		}

		/**
		\brief Facilitates lookup of required endgame type based on tracker type
		
		Your current choices are PSEG or Cauchy.
	
		To get the Power Series Endgame for Adaptive Precision Tracker, use the following example code:
		\code
		using EGT = EndgameSelector<AMPTracker>::PSEG
		\endcode

		\tparam TrackerT The type of tracker you want to use.
		*/
		template<typename TrackerT>
		struct EndgameSelector
		{ };

		template<>
		struct EndgameSelector<DoublePrecisionTracker>
		{
			using PSEG = endgame::FixedPrecPowerSeriesEndgame<DoublePrecisionTracker>;
			using Cauchy = endgame::FixedPrecCauchyEndgame<DoublePrecisionTracker>;
		};

		template<>
		struct EndgameSelector<MultiplePrecisionTracker>
		{
			using PSEG = endgame::FixedPrecPowerSeriesEndgame<MultiplePrecisionTracker>;
			using Cauchy = endgame::FixedPrecCauchyEndgame<MultiplePrecisionTracker>;
		};

		template<class D>
		struct EndgameSelector<FixedPrecisionTracker<D> >
		{
			using PSEG = typename EndgameSelector<D>::PSEG;
			using Cauchy = typename EndgameSelector<D>::Cauchy;
		};

		template<>
		struct EndgameSelector<AMPTracker>
		{
			using PSEG = endgame::AMPPowerSeriesEndgame;
			using Cauchy = endgame::AMPCauchyEndgame;
		};




		enum class SuccessCode
		{
			Success,
			HigherPrecisionNecessary,
			ReduceStepSize,
			GoingToInfinity,
			FailedToConverge,
			MatrixSolveFailure,
			MatrixSolveFailureFirstPartOfPrediction,
			MaxNumStepsTaken,
			MaxPrecisionReached,
			MinStepSizeReached,
			Failure,
			SingularStartPoint,
			ExternallyTerminated,
			MinTrackTimeReached,
			SecurityMaxNormReached,
			CycleNumTooHigh,

		};

		

		/** 
		\namespace config
		*/
		namespace config{

			enum class Predictor
			{
				Constant,
				Euler,
				Heun,
				RK4,
				HeunEuler,
				RKNorsett34,
				RKF45,
				RKCashKarp45,
				RKDormandPrince56,
				RKVerner67
			};

			


			template<typename T>
			struct Tolerances
			{	
				T newton_before_endgame = T(1)/T(100000);
				T newton_during_endgame = T(1)/T(1000000);

				T final_tolerance = T(1)/T(100000000000);
				T final_tolerance_multiplier = T(10); // This multiplier is used to cluster or de-cluster points at the target system. 

				T path_truncation_threshold = T(100000);
				T final_tolerance_times_final_tolerance_multiplier = final_tolerance * final_tolerance_multiplier;
			};


			template<typename T>
			struct Stepping
			{
				T initial_step_size = T(1)/T(10);
				T max_step_size = T(1)/T(10);
				T min_step_size = T(1e-100);

				T step_size_success_factor = T(2);
				T step_size_fail_factor = T(1)/T(2);

				unsigned consecutive_successful_steps_before_stepsize_increase = 5;

				unsigned min_num_steps = 1;
				unsigned max_num_steps = 1e5;

				unsigned frequency_of_CN_estimation = 1;
			};


			
			struct Newton
			{
				unsigned max_num_newton_iterations = 2;
				unsigned min_num_newton_iterations = 1;
			};


			
			struct RenameMe
			{
				
				unsigned sharpendigits;

				mpfr_float function_residual_tolerance;
				mpfr_float ratio_tolerance;
			};

			template<typename T>
			struct Security
			{
				int level = 0;
				T max_norm = T(100000);
			};

			template<typename T>
			struct Endgame
			{
				unsigned num_sample_points = 3;
				T min_track_time = T(1e-100); //nbrh radius in Bertini book.
				T sample_factor = T(1)/T(2);
				unsigned max_num_newton_iterations = 15; // the maximum number allowable iterations during endgames, for points used to approximate the final solution.
			};


			struct PowerSeries
			{
				unsigned max_cycle_number = 6;
				unsigned cycle_number_amplification = 5;
			};

			template<typename T>
			struct Cauchy
			{
				T cycle_cutoff_time = T(1)/T(100000000);
				T ratio_cutoff_time = T(1)/T(100000000000000);
				T minimum_for_c_over_k_stabilization = T(3)/T(4);
				unsigned int num_needed_for_stabilization = 3;
				T maximum_cauchy_ratio = T(1)/T(2);
				unsigned int fail_safe_maximum_cycle_number = 250; //max number of loops before giving up. 

			};


			struct TrackBack
			{
				unsigned minimum_cycle;
				bool junk_removal_test;
				unsigned max_depth_LDT;
			};







			template<typename T>
			struct PostProcessing{
				T real_threshold;
				T endpoint_finite_threshold;
				T final_tol_multiplier;
				T final_tol_times_mult;
			};









			template<typename T>
			struct Regeneration
			{
				bool remove_infinite_endpoints;
				bool higher_dimension_check;

				T newton_before_endgame;
				T newton_during_endgame;
				T final_tolerance;
			};

			




			template<typename ComplexType>
			struct FixedPrecisionConfig
			{
				/**
				\brief Construct a ready-to-go set of fixed precision settings from a system.
				*/
				FixedPrecisionConfig(System const& sys) 
				{ }
			};

			template<typename ComplexType>
			inline
			std::ostream& operator<<(std::ostream & out, FixedPrecisionConfig<ComplexType> const& fpc)
			{
				return out;
			}


			/**
			Holds the program parameters with respect to Adaptive Multiple Precision.
			
			These criteria are developed in \cite AMP1, \cite AMP2.

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

				mpfr_float epsilon;  ///< Bound on growth in error from linear solves.  This is \f$\epsilon\f$ in \cite AMP1, \cite AMP2, and is used for AMP criteria A and B.  See top of page 13 of \cite AMP1.  A pessimistic bound is \f$2^n\f$.
				// rename to linear_solve_error_bound.

				mpfr_float Phi;  ///< Bound on \f$\Phi\f$ (an error bound).   Used for AMP criteria A, B.
				// \f$\Phi\f$ is error in Jacobian evaluation divided by the unit roundoff error, \f$10^{-P}\f$
				// rename to jacobian_eval_error_bound

				mpfr_float Psi;  ///< Bound on \f$\Psi\f$ (an error bound).   Used for AMP criterion C.
				// Error in function evaluation, divided by the precision-dependent unit roundoff error.
				// rename to function_eval_error_bound

				int safety_digits_1 = 1; ///< User-chosen setting for the number of safety digits used during Criteria A & B.
				int safety_digits_2 = 1; ///< User-chosen setting for the number of safety digits used during Criterion C.
				unsigned int maximum_precision = 300; ///< User-chosed setting for the maximum allowable precision.  Paths will die if their precision is requested to be set higher than this threshold.
				
				unsigned consecutive_successful_steps_before_precision_decrease = 10;

				unsigned max_num_precision_decreases = 10; ///< The maximum number of times precision can be lowered during tracking of a segment of path.
				

				/**
				 \brief Set epsilon, degree bound, and coefficient bound from system.
				 
				 * Epsilon is set as the square of the number of variables.
				 * Bound on degree is set from a call to System class.  Let this be \f$D\f$  \see System::DegreeBound().
				 * Bound on absolute values of coeffs is set from a call to System class.  Let this be \f$B\f$.  \see System::CoefficientBound().
				*/
				void SetBoundsAndEpsilonFrom(System const& sys)
				{
					epsilon = pow(mpfr_float(sys.NumVariables()),2);
					degree_bound = sys.DegreeBound();
					coefficient_bound = sys.CoefficientBound();
				}
				

				/**
				 Sets values epsilon, Phi, Psi, degree_bound, and coefficient_bound from input system.
				
				 * Phi becomes \f$ D*(D-1)*B \f$.
				 * Psi is set as \f$ D*B \f$.
				*/
				void SetPhiPsiFromBounds()
				{
					Phi = degree_bound*(degree_bound-mpfr_float(1))*coefficient_bound;
				    Psi = degree_bound*coefficient_bound;  //Psi from the AMP paper.
				}

				void SetAMPConfigFrom(System const& sys)
				{
					SetBoundsAndEpsilonFrom(sys);
					SetPhiPsiFromBounds();
				}

				AdaptiveMultiplePrecisionConfig() : coefficient_bound(1000), degree_bound(5), safety_digits_1(1), safety_digits_2(1), maximum_precision(300) 
				{}

				AdaptiveMultiplePrecisionConfig(System const& sys) : AdaptiveMultiplePrecisionConfig()
				{
					SetAMPConfigFrom(sys);
				}
			}; // re: AdaptiveMultiplePrecisionConfig

			inline
			std::ostream& operator<<(std::ostream & out, AdaptiveMultiplePrecisionConfig const& AMP)
			{
				out << "coefficient_bound: " << AMP.coefficient_bound << "\n";
				out << "degree_bound: " << AMP.degree_bound << "\n";
				out << "epsilon: " << AMP.epsilon << "\n";
				out << "Phi: " << AMP.Phi << "\n";
				out << "Psi: " << AMP.Psi << "\n";
				out << "safety_digits_1: " << AMP.safety_digits_1 << "\n";
				out << "safety_digits_2: " << AMP.safety_digits_2 << "\n";
				out << "consecutive_successful_steps_before_precision_decrease" << AMP.consecutive_successful_steps_before_precision_decrease << "\n";
				return out;
			}

			
			/**
			\brief Construct a ready-to-go set of AMP settings from a system.
			
			

			\see AdaptiveMultiplePrecisionConfig::SetBoundsAndEpsilonFrom
			\see AdaptiveMultiplePrecisionConfig::SetPhiPsiFromBounds
			\see AdaptiveMultiplePrecisionConfig::SetAMPConfigFrom
			*/
			inline
			static
			AdaptiveMultiplePrecisionConfig AMPConfigFrom(System const& sys) 
			{
				AdaptiveMultiplePrecisionConfig AMP;				
				AMP.SetAMPConfigFrom(sys);
			    return AMP;
			}

		} //re: namespace config
	} // re: namespace tracking 
} // re: namespace bertini


#endif
