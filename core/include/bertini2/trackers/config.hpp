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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

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


    /// \brief The precision regime a tracker operates in: fixed double, fixed multiple, or adaptive.
    enum class PrecisionType //E.2.1
    {
        Fixed,
        FixedMultiple,
        Adaptive
    };


    /// \brief The predictor (ODE integration) method used during path tracking.
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







    /**
    \brief Metadata produced by a single predict or correct step.

    Collapses the formerly hand-threaded out-parameters (norm_J, norm_J_inverse,
    condition_number_estimate, size_proportion, error_estimate, norm_delta_z) into one
    struct.  Not every field is written by every step: size_proportion/error_estimate are
    predictor-only (and error_estimate only for embedded methods); norm_delta_z is
    corrector-only.  Unwritten fields keep their default of 0.
    */
    struct StepMetadata
    {
        NumErrorT norm_J = 0;                    ///< ||J|| (Frobenius) at the step.
        NumErrorT norm_J_inverse = 0;            ///< estimate of ||J^{-1}|| via the condition probe.
        NumErrorT condition_number_estimate = 0; ///< norm_J * norm_J_inverse (refreshed per frequency_of_CN_estimation).
        NumErrorT size_proportion = 0;           ///< AMP "a" (predictor only).
        NumErrorT error_estimate = 0;            ///< embedded-method error estimate (predictor, embedded only).
        NumErrorT norm_delta_z = 0;              ///< ||latest Newton step|| (corrector only).
    };


    /// \brief Settings governing step-size adjustment during tracking.
    struct SteppingConfig
    {
        // mpq_rational: exact rationals with no MPFR precision state — safe in DefaultConstruct<T>::value statics.
        // real_mp fields here would be initialized at BMP's startup precision (20) and contaminate
        // tracker arithmetic when target precision < 20 via preserve_related_precision.
        mpq_rational initial_step_size{1, 10}; ///< The length of the first time step when calling TrackPath.  StepInitSize
        mpq_rational max_step_size{1, 10};     ///<  The largest allowed step size.  MaxStepSize
        double       min_step_size = 1e-100;   ///< The minimum allowed step size (threshold only, double precision is sufficient).  MinStepSize

        mpq_rational step_size_success_factor{2, 1}; ///< Factor by which to dilate the time step when triggered.  StepSuccessFactor
        mpq_rational step_size_fail_factor{1, 2};    ///< Factor by which to contract the time step when triggered.  StepFailFactor

        unsigned consecutive_successful_steps_before_stepsize_increase = 5; ///< What it says.  If you can come up with a better name, please suggest it.  StepsForIncrease

        unsigned min_num_steps = 1; ///< The minimum number of steps allowed during tracking.
        unsigned max_num_steps = 1e5; ///< The maximum number of steps allowed during tracking, counting FAILED steps as well as successful ones.  This is per call to TrackPath.  MaxNumberSteps

        unsigned frequency_of_CN_estimation = 1; ///< Estimate the condition number every so many steps.  Eh.
    };



    /// \brief Settings governing the Newton corrector's iteration bounds.
    struct NewtonConfig
    {
        unsigned max_num_newton_iterations = 2; ///< The maximum number of Newton iterations per correction.  MaxNewtonIts
        unsigned min_num_newton_iterations = 1; ///< The minimum number of Newton iterations per correction.
    };








    /// \brief Settings for a fixed-precision tracker (carries the single working precision in effect).
    struct FixedPrecisionConfig
    {
        using RealT = double;  ///< The real number type.

        /**
        \brief The number of digits to always work at.

        For a double-precision tracker this is DoublePrecision() (16) and cannot be changed.  For a
        fixed-multiple tracker it is the precision the whole solve runs at -- the tracker, the system,
        the start points, and the working precision all sit at this one value.  A tracker keeps this
        field in sync with its actual precision, so reading it tells you the precision in effect; set it
        to choose a different fixed precision (the algorithm then lifts the system and start points to
        match).  The sentinel 0 means "unset -- use the tracker's natural precision".
        */
        unsigned precision = 0;

        /**
        \brief Construct a ready-to-go set of fixed precision settings from a system.
        */
        explicit
        FixedPrecisionConfig(System const& /*sys*/)
        { }

        FixedPrecisionConfig() = default;
    };


    /// \brief Stream-insertion for FixedPrecisionConfig (a no-op; the config carries no printable state).
    inline
    std::ostream& operator<<(std::ostream & out, FixedPrecisionConfig const& /*fpc*/)
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

        /// \brief Bound on the growth in error from a linear solve.  Used for AMP criteria A and
        /// B, where it is called \f$\epsilon\f$ in \cite AMP1 and \cite AMP2; see the top of page
        /// 13 of the former.  A pessimistic bound is \f$2^n\f$; we use \f$n^2\f$ in the number of
        /// variables.
        NumErrorT linear_solve_error_bound;

        /// \brief Bound on the error of evaluating the Jacobian, divided by the unit roundoff
        /// error \f$10^{-P}\f$.  Used for AMP criteria A and B, where it is called \f$\Phi\f$
        /// in \cite AMP1 and \cite AMP2.
        ///
        /// For a polynomial system this is \f$D(D-1)B\f$ in the degree and coefficient bounds
        /// below; that recipe is what SetErrorBoundsFromDegreeAndCoefficient applies.  The
        /// quantity itself is more general than the recipe, which is why it is spelled out
        /// rather than named after a Greek letter: a caller with a system that has no degree
        /// can set this directly.
        NumErrorT jacobian_eval_error_bound;

        /// \brief Bound on the error of evaluating the functions, divided by the unit roundoff
        /// error.  Used for AMP criterion C, where it is called \f$\Psi\f$ in \cite AMP1 and
        /// \cite AMP2.
        ///
        /// For a polynomial system this is \f$DB\f$; see the note above.
        NumErrorT function_eval_error_bound;

        int safety_digits_1 = 1; ///< User-chosen setting for the number of safety digits used during Criteria A & B.
        int safety_digits_2 = 1; ///< User-chosen setting for the number of safety digits used during Criterion C.
        unsigned int maximum_precision = 300; ///< User-chosed setting for the maximum allowable precision.  Paths will die if their precision is requested to be set higher than this threshold.

        // Note: a single setting -- Bertini 1's StepsForIncrease, i.e.
        // SteppingConfig::consecutive_successful_steps_before_stepsize_increase -- gates BOTH stepsize
        // increase AND precision decrease (the required number of consecutive successful steps).
        // Precision decrease is additionally subject to B1's digits-margin hysteresis (see
        // ExtraDigitsBeforePrecisionDecrease in amp_tracker.hpp).  The old, duplicate AMP-config setting
        // `consecutive_successful_steps_before_precision_decrease` (a B2 deviation) was removed.

        unsigned max_num_precision_decreases = 10; ///< The maximum number of times precision can be lowered during tracking of a segment of path.


        /**
         \brief Set the linear-solve error bound, the degree bound and the coefficient bound from
         a system.

         * The linear-solve error bound is the square of the number of variables.
         * Bound on degree is set from a call to System class.  Let this be \f$D\f$  \see System::DegreeBound().
         * Bound on absolute values of coeffs is set from a call to System class.  Let this be \f$B\f$.  \see System::CoefficientBound().
        */
        void SetBoundsFrom(System const& sys)
        {
            using std::pow;

            // Refuse before DegreeBound() does, so the message is about what the caller was
            // trying to do rather than about degrees.  The two bounds are on the error of
            // evaluating the Jacobian and the functions; the degree and coefficient bounds are
            // merely the recipe for them that holds for polynomials.  There is no accepted
            // recipe for an analytic system -- see issue #439 -- so we do not invent one.
            if (!sys.IsPolynomial())
                throw std::runtime_error("adaptive precision derives its error bounds from the "
                    "degree of the system, and this system is not a polynomial one, so it has no "
                    "degree.  Two ways on: track with a fixed-precision tracker, at double or at "
                    "multiple precision, which needs no such bound; or set this config's "
                    "jacobian_eval_error_bound and function_eval_error_bound yourself and hand it "
                    "to the tracker, choosing values you can defend for your system.");

            linear_solve_error_bound = pow(NumErrorT(sys.NumVariables()),2);
            degree_bound = sys.DegreeBound();
            coefficient_bound = sys.CoefficientBound<complex_dbl>();
        }


        /// \brief Set the two evaluation error bounds from the degree and coefficient bounds.
        ///
        /// The polynomial recipe from \cite AMP1: with \f$D\f$ the degree bound and \f$B\f$ the
        /// coefficient bound, the Jacobian evaluation error is bounded by \f$D(D-1)B\f$ and the
        /// function evaluation error by \f$DB\f$.  It applies only to polynomial systems, which
        /// is the only kind that has a \f$D\f$; a caller with an analytic system sets the two
        /// bounds directly instead.
        void SetErrorBoundsFromDegreeAndCoefficient()
        {
            jacobian_eval_error_bound = degree_bound*(degree_bound-NumErrorT(1))*coefficient_bound;
            function_eval_error_bound = degree_bound*coefficient_bound;
        }

        /// \brief Set every system-derived AMP setting from a system: the degree, coefficient and
        /// linear-solve bounds, and the two evaluation error bounds derived from them.
        void SetAMPConfigFrom(System const& sys)
        {
            SetBoundsFrom(sys);
            SetErrorBoundsFromDegreeAndCoefficient();
        }

        /// \brief Construct with default AMP bounds and safety digits.
        ///
        /// The three error bounds are normally recomputed from the system before tracking
        /// (\see SetAMPConfigFrom).  They are nonetheless initialized here so that a
        /// default-constructed config is fully deterministic: the linear-solve bound takes the
        /// single-variable value (\f$1^2\f$), and the two evaluation bounds are derived from the
        /// default degree and coefficient bounds so the object is internally consistent.  Leaving
        /// them uninitialized produced garbage (e.g. NaN) that broke value comparison and pickle
        /// round-tripping.
        AdaptiveMultiplePrecisionConfig() : coefficient_bound(1000), degree_bound(5), linear_solve_error_bound(1), safety_digits_1(1), safety_digits_2(1), maximum_precision(300)
        {
            SetErrorBoundsFromDegreeAndCoefficient();
        }

        /// \brief Construct AMP settings derived from a system's bounds, where that is possible.
        ///
        /// A system that is not polynomial has no degree bound, so there is nothing to derive from
        /// and the defaults are kept.  This constructor is what a solver's automatic setup calls,
        /// and it must not refuse, or an analytic homotopy could not be handed a configuration the
        /// caller chose: the refusal would fire first, during construction.  The defaults are
        /// internally consistent and are NOT a claim about such a system, so
        /// \see SetAMPConfigFrom, which a deliberate caller uses, does refuse.  A caller who wants
        /// adaptive precision here supplies the two evaluation error bounds themselves.
        explicit
        AdaptiveMultiplePrecisionConfig(System const& sys) : AdaptiveMultiplePrecisionConfig()
        {
            if (sys.IsPolynomial())
                SetAMPConfigFrom(sys);
        }
    }; // re: AdaptiveMultiplePrecisionConfig

    /// \brief Stream-insertion for AdaptiveMultiplePrecisionConfig, printing its bounds and safety digits.
    inline
    std::ostream& operator<<(std::ostream & out, AdaptiveMultiplePrecisionConfig const& AMP)
    {
        out << "coefficient_bound: " << AMP.coefficient_bound << "\n";
        out << "degree_bound: " << AMP.degree_bound << "\n";
        out << "linear_solve_error_bound: " << AMP.linear_solve_error_bound << "\n";
        out << "jacobian_eval_error_bound: " << AMP.jacobian_eval_error_bound << "\n";
        out << "function_eval_error_bound: " << AMP.function_eval_error_bound << "\n";
        out << "safety_digits_1: " << AMP.safety_digits_1 << "\n";
        out << "safety_digits_2: " << AMP.safety_digits_2 << "\n";
        out << "max_num_precision_decreases: " << AMP.max_num_precision_decreases << "\n";
        return out;
    }


    /**
    \brief Construct a ready-to-go set of AMP settings from a system.



    \see AdaptiveMultiplePrecisionConfig::SetBoundsFrom
    \see AdaptiveMultiplePrecisionConfig::SetErrorBoundsFromDegreeAndCoefficient
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
    /// \brief Trait lookup mapping a tracker type to its numeric types, event-emitter type, precision
    ///        config, and the type/config lists it needs.  Specialized per concrete tracker type.
    template<class T>
    struct TrackerTraits
    {};


    /// \cond TRACKER_TRAITS_SPECIALIZATIONS

    template<>
    struct TrackerTraits<DoublePrecisionTracker>
    {
        static constexpr char const* kRecordName = "double";  ///< Stable tracker name for records (b2rec ask identity; never typeid).
        using BaseComplexT = complex_dbl;
        using BaseRealT = double;
        using EventEmitterType = FixedPrecisionTracker<DoublePrecisionTracker>;
        using PrecisionConfig = FixedPrecisionConfig;
        enum {
            IsFixedPrec = 1,
            IsAdaptivePrec = 0
        };

        using NeededTypes = detail::TypeList<complex_dbl>;
        using NeededConfigs = detail::TypeList<
            SteppingConfig,
            NewtonConfig,
            PrecisionConfig
            >;
    };


    template<>
    struct TrackerTraits<MultiplePrecisionTracker>
    {
        static constexpr char const* kRecordName = "multiple";  ///< Stable tracker name for records (b2rec ask identity; never typeid).
        using BaseComplexT = complex_mp;
        using BaseRealT = real_mp;
        using EventEmitterType = FixedPrecisionTracker<MultiplePrecisionTracker>;
        using PrecisionConfig = FixedPrecisionConfig;

        enum {
            IsFixedPrec = 1,
            IsAdaptivePrec = 0
        };

        using NeededTypes = detail::TypeList<complex_mp>;

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
        static constexpr char const* kRecordName = "adaptive";  ///< Stable tracker name for records (b2rec ask identity; never typeid).
        using BaseComplexT = complex_mp;
        using BaseRealT = real_mp;
        using EventEmitterType = AMPTracker;
        using PrecisionConfig = AdaptiveMultiplePrecisionConfig;

        enum {
            IsFixedPrec = 0,
            IsAdaptivePrec = 1
        };

        using NeededTypes = detail::TypeList<complex_dbl, complex_mp>;

        using NeededConfigs = detail::TypeList<
            SteppingConfig,
            NewtonConfig,
            PrecisionConfig
            >;
    };




    template<class D>
    struct TrackerTraits<FixedPrecisionTracker<D> > : public TrackerTraits<D>
    {
        using BaseComplexT = typename TrackerTraits<D>::BaseComplexT;
        using BaseRealT = typename TrackerTraits<D>::BaseRealT;
        using EventEmitterType = typename TrackerTraits<D>::EventEmitterType;
        using PrecisionConfig = typename TrackerTraits<D>::PrecisionConfig;

        enum {
            IsFixedPrec = 0,
            IsAdaptivePrec = 1
        };

        using NeededTypes = typename TrackerTraits<D>::NeededTypes;
        using NeededConfigs = typename TrackerTraits<D>::NeededConfigs;
    };

    /// \endcond

} // re: namespace tracking
} // re: namespace bertini


#endif
