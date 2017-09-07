//This file is part of Bertini 2.
//
//trackers/include/bertini2/trackers/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//trackers/include/bertini2/trackers/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/include/bertini2/trackers/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University

#ifndef BERTINI_TRACKING_CONFIG_HPP
#define BERTINI_TRACKING_CONFIG_HPP

/**
\file include/bertini2/trackers/config.hpp

\brief Configs and settings for tracking.
*/
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/detail/typelist.hpp"

#include "bertini2/common/config.hpp"


namespace bertini
{

namespace tracking{

	
	enum class PrecisionType //E.2.1
	{
		Fixed,
		Adaptive
	};
	

	enum class Predictor //E.4.3
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

	


	


	struct SteppingConfig
	{
		using T = mpq_rational;

		T initial_step_size = T(1)/T(10); ///< The length of the first time step when calling TrackPath.  You can turn it resetting, so subsequent calls use the same stepsize, too.  You make a call to the Tracker itself.
		T max_step_size = T(1)/T(10); ///<  The largest allowed step size.  MaxStepSize
		T min_step_size = T(1)/T(1e100); ///< The mimum allowed step size.  MinStepSize

		T step_size_success_factor = T(2); ///< Factor by which to dilate the time step when triggered.  StepSuccessFactor
		T step_size_fail_factor = T(1)/T(2); ///< Factor by which to contract the time step when triggered.  StepFailFactor

		unsigned consecutive_successful_steps_before_stepsize_increase = 5; ///< What it says.  If you can come up with a better name, please suggest it.  StepsForIncrease

		unsigned min_num_steps = 1; ///< The minimum number of steps allowed during tracking.
		unsigned max_num_steps = 1e5; ///< The maximum number of steps allowed during tracking.  This is per call to TrackPath.  MaxNumberSteps

		unsigned frequency_of_CN_estimation = 1; ///< Estimate the condition number every so many steps.  Eh.
	};


	
	struct NewtonConfig
	{
		unsigned max_num_newton_iterations = 2; //MaxNewtonIts
		unsigned min_num_newton_iterations = 1;
	};


	
	

	


	struct FixedPrecisionConfig
	{
		using RealType = double;

		/**
		\brief Construct a ready-to-go set of fixed precision settings from a system.
		*/
		explicit
		FixedPrecisionConfig(System const& sys) 
		{ }

		FixedPrecisionConfig() = default;
	};


	inline
	std::ostream& operator<<(std::ostream & out, FixedPrecisionConfig const& fpc)
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
		NumErrorT coefficient_bound;  ///< User-defined bound on the sum of the abs vals of the coeffs for any polynomial in the system (for adaptive precision). 
		NumErrorT degree_bound; ///<  User-set bound on degrees of polynomials in the system - tricky to compute for factored polys, subfuncs, etc. (for adaptive precision). 

		NumErrorT epsilon;  ///< Bound on growth in error from linear solves.  This is \f$\epsilon\f$ in \cite AMP1, \cite AMP2, and is used for AMP criteria A and B.  See top of page 13 of \cite AMP1.  A pessimistic bound is \f$2^n\f$.
		// rename to linear_solve_error_bound.

		NumErrorT Phi;  ///< Bound on \f$\Phi\f$ (an error bound).   Used for AMP criteria A, B.
		// \f$\Phi\f$ is error in Jacobian evaluation divided by the unit roundoff error, \f$10^{-P}\f$
		// rename to jacobian_eval_error_bound

		NumErrorT Psi;  ///< Bound on \f$\Psi\f$ (an error bound).   Used for AMP criterion C.
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
			using std::pow;

			epsilon = pow(NumErrorT(sys.NumVariables()),2);
			degree_bound = sys.DegreeBound();
			coefficient_bound = sys.CoefficientBound<dbl>();
		}
		

		/**
		 Sets values epsilon, Phi, Psi, degree_bound, and coefficient_bound from input system.
		
		 * Phi becomes \f$ D*(D-1)*B \f$.
		 * Psi is set as \f$ D*B \f$.
		*/
		void SetPhiPsiFromBounds()
		{	
			Phi = degree_bound*(degree_bound-NumErrorT(1))*coefficient_bound;
		    Psi = degree_bound*coefficient_bound;  //Psi from the AMP paper.
		}

		void SetAMPConfigFrom(System const& sys)
		{
			SetBoundsAndEpsilonFrom(sys);
			SetPhiPsiFromBounds();
		}

		AdaptiveMultiplePrecisionConfig() : coefficient_bound(1000), degree_bound(5), safety_digits_1(1), safety_digits_2(1), maximum_precision(300) 
		{}

		explicit
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

	// forward declarations
	template<class D>
	class Tracker;
	template<class D>
	class FixedPrecisionTracker;
	class MultiplePrecisionTracker;
	class DoublePrecisionTracker;
	

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
		using PrecisionConfig = FixedPrecisionConfig;
		enum {
			IsFixedPrec = 1,
			IsAdaptivePrec = 0
		};

		using NeededTypes = detail::TypeList<dbl>;
		using NeededConfigs = detail::TypeList<
			SteppingConfig, 
			NewtonConfig,
			PrecisionConfig
			>;
	};


	template<>
	struct TrackerTraits<MultiplePrecisionTracker>
	{
		using BaseComplexType = mpfr;
		using BaseRealType = mpfr_float;
		using EventEmitterType = FixedPrecisionTracker<MultiplePrecisionTracker>;
		using PrecisionConfig = FixedPrecisionConfig;

		enum {
			IsFixedPrec = 1,
			IsAdaptivePrec = 0
		};

		using NeededTypes = detail::TypeList<mpfr>;

		using NeededConfigs = detail::TypeList<
			SteppingConfig, 
			NewtonConfig,
			PrecisionConfig
			>;
	};


	class AMPTracker; // forward declare
	template<>
	struct TrackerTraits<AMPTracker>
	{
		using BaseComplexType = mpfr;
		using BaseRealType = mpfr_float;
		using EventEmitterType = AMPTracker;
		using PrecisionConfig = AdaptiveMultiplePrecisionConfig;

		enum {
			IsFixedPrec = 0,
			IsAdaptivePrec = 1
		};

		using NeededTypes = detail::TypeList<dbl, mpfr>;

		using NeededConfigs = detail::TypeList<
			SteppingConfig, 
			NewtonConfig,
			PrecisionConfig
			>;
	};


	

	template<class D>
	struct TrackerTraits<FixedPrecisionTracker<D> > : public TrackerTraits<D>
	{ 
		using BaseComplexType = typename TrackerTraits<D>::BaseComplexType;
		using BaseRealType = typename TrackerTraits<D>::BaseRealType;
		using EventEmitterType = typename TrackerTraits<D>::EventEmitterType;
		using PrecisionConfig = typename TrackerTraits<D>::PrecisionConfig;

		enum {
			IsFixedPrec = 0,
			IsAdaptivePrec = 1
		};

		using NeededTypes = typename TrackerTraits<D>::NeededTypes;
		using NeededConfigs = typename TrackerTraits<D>::NeededConfigs;
	};

} // re: namespace tracking 
} // re: namespace bertini


#endif
