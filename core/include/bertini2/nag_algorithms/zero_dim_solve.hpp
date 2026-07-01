//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/zero_dim_solve.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/zero_dim_solve.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/nag_algorithms/zero_dim_solve.hpp

\brief Provides the algorithm for computing all zero-dimensional solutions for an algberaic system.
*/

#pragma once

#include "bertini2/num_traits.hpp"

#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"
#include "bertini2/nag_algorithms/midpath_check.hpp"
#include "bertini2/io/generators.hpp"

#include "bertini2/detail/configured.hpp"
#include "bertini2/detail/observable.hpp"

#include "bertini2/nag_algorithms/common/algorithm_base.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"
#include "bertini2/nag_algorithms/events.hpp"
#include "bertini2/system/start_base.hpp"   // start_system::StartSystem + StartSystemFactory / MakeStartFactory
#include "bertini2/parallel.hpp"
#include <chrono>
#include <mutex>
#include <iostream>
#include <map>


namespace bertini {

	// forward-declare the interim default start system so ZeroDimSolver's default factory argument
	// (MakeStartFactory<TotalDegreeBinomial>) can name it; the concrete type rides in via start_systems.hpp
	// at every call site that actually constructs a ZeroDimSolver.
	namespace start_system { class TotalDegreeBinomial; }

	namespace algorithm {


/**
\brief The continuation primitive: given a homotopy + a source of start points, track each path
through the tracker and endgame, resolve crossings, classify the endpoints, and report.

This is the engine.  It holds the homotopy, the start system (which supplies the start points), and
a target system (for dehomogenize / residual / classification) by *reference* -- the caller owns
them.  ZeroDimSolver is the algorithm built on top: it owns and builds those systems, then drives
this engine (see below).
*/
template<typename TrackerType, typename EndgameType, typename SystemType>
struct HomotopySolver;

/**
\brief The zero-dimensional solve algorithm: given a polynomial system, form a start system and a
homotopy, then run the continuation engine.  Owns its systems; is-a HomotopySolver.
*/
template<typename TrackerType, typename EndgameType, typename SystemType>
struct ZeroDimSolver;



/**
specify the traits for the algorithm.  this is why we need the forward declare.  ZeroDimSolver
derives from HomotopySolver (and so inherits its Configured base), so only HomotopySolver needs
traits specialized here.
*/
template<typename TrackerType, typename EndgameType, typename SystemType>
struct AlgoTraits <HomotopySolver<TrackerType, EndgameType, SystemType>>
{
	using BaseRealT = typename tracking::TrackerTraits<TrackerType>::BaseRealT;        ///< The real number type of the tracker.
	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;  ///< The complex number type of the tracker.

	/// The config types this algorithm reads (drives the reusable Python config interface).
	using NeededConfigs = detail::TypeList<
								TolerancesConfig,
								PostProcessingConfig,
								ZeroDimConfig,
								AutoRetrackConfig
								>;
};


/// \brief Type-erased base for a zero-dimensional solve, exposing output and config entry points.
struct AnyZeroDim : public virtual AnyAlgorithm
{
	/// \brief Write the human-readable main-data report to \p out.
	virtual void WriteMainData(std::ostream& out) const = 0;
	/// \brief Write the raw-data report to \p out.
	virtual void WriteRawData(std::ostream& out)  const = 0;
	// Bertini 1.7-compatible machine-readable solution files (count-led coordinate blocks).
	/// \brief Write the finite solutions (Bertini 1.7 format) to \p out.
	virtual void WriteFiniteSolutions(std::ostream& out)      const = 0;
	/// \brief Write the real finite solutions (Bertini 1.7 format) to \p out.
	virtual void WriteRealFiniteSolutions(std::ostream& out)  const = 0;
	/// \brief Write the nonsingular solutions (Bertini 1.7 format) to \p out.
	virtual void WriteNonsingularSolutions(std::ostream& out) const = 0;
	/// \brief Write the singular solutions (Bertini 1.7 format) to \p out.
	virtual void WriteSingularSolutions(std::ostream& out)    const = 0;
	/// \brief Write the raw solutions (Bertini 1.7 format) to \p out.
	virtual void WriteRawSolutions(std::ostream& out)         const = 0;
	/// \brief Apply configuration settings parsed from a classic-input config string.
	virtual void ApplyParsedConfigs(std::string const& config_str) = 0;
	virtual ~AnyZeroDim() = default;
};




/// Index type for solutions/paths (over a double-precision solution container).
using SolnIndT = typename SolnCont<complex_dbl>::size_type;


/// metadata structs

/// \brief Solve-wide bookkeeping: path counts and timing.
struct AlgorithmMetaData
{

	SolnIndT number_path_failures = 0;   ///< Number of paths that failed.
	SolnIndT number_path_successes = 0;  ///< Number of paths that succeeded.
	SolnIndT number_paths_tracked = 0;   ///< Total number of paths tracked.

	std::chrono::system_clock::time_point start_time;  ///< Wall-clock time the solve started.
	std::chrono::microseconds elapsed_time;            ///< Total wall-clock time the solve took.
};


/// \brief Per-solution metadata gathered during and after a path track.
template<typename ComplexT>
struct SolutionMetaData
{
	using SolnIndT = typename SolnCont<ComplexT>::size_type;  ///< Index type for solutions/paths.

	// only vaguely metadata.  artifacts of randomness or ordering
	SolnIndT path_index;     		///< Path number of the solution.
	SolnIndT solution_index;      	///< Solution number.

	///// things computed across all of the solve
	bool precision_changed = false;          ///< Whether the working precision changed during the path.
	ComplexT time_of_first_prec_increase;    ///< Time value of the first increase in precision.
	decltype(DefaultPrecision()) max_precision_used = 0;  ///< Highest precision used on the path.
	double path_time_seconds = 0.0;          ///< Wall-clock time for the whole path (pre-endgame + endgame), seconds.  Not an identity field (excluded from operator==).

	///// things computed in pre-endgame only
	SuccessCode pre_endgame_success_code = SuccessCode::NeverStarted;     ///< Success code of the pre-endgame track.


	///// things computed in endgame only
	NumErrorT condition_number; 				///< The latest estimate of the condition number.
	NumErrorT newton_residual; 				///< The latest Newton residual.
	ComplexT final_time_used;   			///< The final value of time tracked to.
	NumErrorT accuracy_estimate; 			///< Accuracy estimate between extrapolations.
	NumErrorT accuracy_estimate_user_coords;	///< Accuracy estimate between extrapolations, in natural coordinates.
	unsigned cycle_num;    						///< Cycle number used in extrapolations.
	// Honest precision/accuracy of the computed solution, as DIGIT COUNTS.  precision_digits is the
	// working precision the endgame actually finished in (DoublePrecision() ~16 for a path that stayed in
	// the adaptive-numeric-type endgame's hardware-double fast lane; the mpfr precision for one that
	// escalated).  accuracy_digits is how many of those digits are trustworthy, from the convergence
	// agreement: floor(-log10(accuracy_estimate)), clamped to [0, precision_digits].  Read together:
	// "computed in precision_digits digits, good to accuracy_digits of them."
	unsigned precision_digits = 0;  ///< Working precision (digits) the endgame finished this solution in.
	unsigned accuracy_digits = 0;   ///< Trustworthy digit count, from the convergence agreement.
	SuccessCode endgame_success_code = SuccessCode::NeverStarted;      ///< Success code of the endgame.


	///// things added by post-processing
	NumErrorT function_residual; 	///< The latest function residual.

	int multiplicity = 1; 		///< Multiplicity of the solution.
	bool multiplicity_representative = true; ///< Whether this is the chosen single representative of its multiplicity cluster (the m copies share one point; exactly one is the representative).
	bool is_real = false;       		///< Whether the (dehomogenized) endpoint is real.
	bool is_finite = false;     		///< Whether the endpoint is finite (not at infinity).
	bool is_singular = false;       		///< Whether the endpoint is singular (multiple, or ill-conditioned).
	// nonsolution flag: a finite, successful endpoint that is NOT a solution of the actual target
	// system -- a nonsolution.  ZeroDimSolver sets this when it squares up an over-determined system:
	// the randomized square system has extraneous roots that satisfy the random combinations but not
	// the original equations.  Orthogonal to is_finite (a nonsolution is finite); the finite / real /
	// singular accessors exclude nonsolutions, and they are exposed on their own (Nonsolutions()).
	// Load-bearing for the regeneration cascade, which must identify and discard nonsolutions.
	bool is_nonsolution = false;  ///< Whether this finite, successful endpoint is not a solution of the actual target system (an extraneous root from squaring up an over-determined system).

	/// \brief Equality comparison over the identity fields (path_time_seconds is excluded).
	bool operator==(const SolutionMetaData<ComplexT> & other){
		bool result =
			this->path_index == other.path_index
			 && this->solution_index == other.solution_index
			 && this->precision_changed == other.precision_changed
			 && this->time_of_first_prec_increase == other.time_of_first_prec_increase
			 && this->max_precision_used == other.max_precision_used
			 && this->pre_endgame_success_code == other.pre_endgame_success_code
			 && this->condition_number == other.condition_number
			 && this->newton_residual == other.newton_residual
			 && this->final_time_used == other.final_time_used
			 && this->accuracy_estimate == other.accuracy_estimate
			 && this->accuracy_estimate_user_coords == other.accuracy_estimate_user_coords
			 && this->cycle_num == other.cycle_num
			 && this->precision_digits == other.precision_digits
			 && this->accuracy_digits == other.accuracy_digits
			 && this->endgame_success_code == other.endgame_success_code
			 && this->function_residual == other.function_residual
			 && this->multiplicity == other.multiplicity
			 && this->multiplicity_representative == other.multiplicity_representative
			 && this->is_real == other.is_real
			 && this->is_finite == other.is_finite
			 && this->is_singular == other.is_singular
			 && this->is_nonsolution == other.is_nonsolution
		;

		return result; }
};

/// \brief Stream insertion: write the metadata fields one per line (used by the Python bindings).
// this is for interoperability with vectors of these in the Python bindings, for better or for worse.
template<typename NumT>
std::ostream& operator<<(std::ostream & out, const SolutionMetaData<NumT> & meta){
	out << "path_index = " << meta.path_index << std::endl;
	out << "solution_index = " << meta.solution_index << std::endl;

	out << "precision_changed = " << meta.precision_changed << std::endl;
	out << "time_of_first_prec_increase = " << meta.time_of_first_prec_increase << std::endl;
	out << "max_precision_used = " << meta.max_precision_used << std::endl;
	out << "path_time_seconds = " << meta.path_time_seconds << std::endl;

	out << "pre_endgame_success_code = " << meta.pre_endgame_success_code << std::endl;

	out << "condition_number = " << meta.condition_number << std::endl;
	out << "newton_residual = " << meta.newton_residual << std::endl;
	out << "final_time_used = " << meta.final_time_used << std::endl;
	out << "accuracy_estimate = " << meta.accuracy_estimate << std::endl;
	out << "accuracy_estimate_user_coords = " << meta.accuracy_estimate_user_coords << std::endl;
	out << "precision_digits = " << meta.precision_digits << std::endl;
	out << "accuracy_digits = " << meta.accuracy_digits << std::endl;
	out << "cycle_num = " << meta.cycle_num << std::endl;
	out << "endgame_success_code = " << meta.endgame_success_code << std::endl;

	out << "function_residual = " << meta.function_residual << std::endl;

	out << "multiplicity = " << meta.multiplicity << std::endl;
	out << "multiplicity_representative = " << meta.multiplicity_representative << std::endl;
	out << "is_real = " << meta.is_real << std::endl;
	out << "is_finite = " << meta.is_finite << std::endl;
	out << "is_singular = " << meta.is_singular << std::endl;
	out << "is_nonsolution = " << meta.is_nonsolution << std::endl;

	return out;
}


/// \brief State captured at the endgame boundary, handed off to start the endgame.
template<typename ComplexT>
struct EGBoundaryMetaData
{
	using RealT = typename NumTraits<ComplexT>::Real;  ///< The real number type.

	Vec<ComplexT> path_point;                            ///< The path point at the endgame boundary.
	SuccessCode success_code = SuccessCode::NeverStarted; ///< Success code of the pre-endgame track.
	RealT last_used_stepsize;                            ///< The last step size used before the boundary.
	// The precision the tracker was actually using when it reached the endgame boundary.
	// Carried explicitly (rather than inferred from Precision(path_point)) because the
	// path_point is always stored as the tracker's BaseComplexT (multiprecision for AMP);
	// a path that tracked in double gets widened on output, so its mantissa precision no
	// longer reflects the precision that was in use.  The endgame should resume at this
	// precision.  See zero_dim_solve TrackSinglePathDuringEG.
	unsigned precision = DoublePrecision();  ///< The precision the tracker was using at the endgame boundary (carried explicitly; see note above).

	/// \cond
	EGBoundaryMetaData() = default;
	EGBoundaryMetaData(EGBoundaryMetaData const&) = default;
	EGBoundaryMetaData& operator=(EGBoundaryMetaData const&) = default;
	/// \endcond
	/// \brief Construct from the boundary point, its success code, the last step size, and the precision in use.
	EGBoundaryMetaData(Vec<ComplexT> const& pt, SuccessCode const& code, RealT const& ss, unsigned prec) :
		path_point(pt), success_code(code), last_used_stepsize(ss), precision(prec)
	{}

	/// \brief Equality comparison (compares all fields).
	bool operator==(const EGBoundaryMetaData<ComplexT> & other){
		bool result =
			this->path_point == other.path_point
			&& this->success_code == other.success_code
			&& this->last_used_stepsize == other.last_used_stepsize
			&& this->precision == other.precision
		;

		return result;
	}
};

/// \brief Stream insertion: write the boundary metadata fields one per line (used by the Python bindings).
// this is for interoperability with vectors of these in the Python bindings, for better or for worse.
template<typename NumT>
std::ostream& operator<<(std::ostream & out, const EGBoundaryMetaData<NumT> & meta){
	out << "path_point = " << meta.path_point << std::endl;
	out << "success_code = " << meta.success_code << std::endl;
	out << "last_used_stepsize = " << meta.last_used_stepsize << std::endl;
	out << "precision = " << meta.precision << std::endl;
	return out;
}

/**
\brief Summary of what the midpath (path-crossing) check found at the endgame boundary, and what the
algorithm did about it.

Two distinct paths landing on the same point at the endgame boundary is a probability-0 event: it
signals a path crossing (under-resolved tracking), not a benign coincidence.  The zero-dim algorithm
detects these at the boundary using a relaxed same-point tolerance and re-tracks the offending paths
with tightened settings (see ZeroDim::EGBoundaryAction).  This struct records the outcome so a caller
can tell whether a solve hit crossings and whether they were resolved.

\see ZeroDim::EndgameBoundaryMetadata
*/
struct MidpathCheckReport
{
	bool passed = true;                                    ///< Did the final midpath check pass (no crossings remained)?
	unsigned num_crossings_detected = 0;                   ///< Number of crossed paths found on the *first* check, before any re-tracking.
	unsigned num_resolve_attempts = 0;                     ///< How many re-track attempts were actually performed.
	std::vector<unsigned long long> crossed_path_indices;  ///< Indices of the paths flagged as crossed on the first check.
};

/// \brief Stream insertion: write the midpath-check report fields one per line.
inline
std::ostream& operator<<(std::ostream & out, const MidpathCheckReport & r){
	out << "passed = " << std::boolalpha << r.passed << std::endl;
	out << "num_crossings_detected = " << r.num_crossings_detected << std::endl;
	out << "num_resolve_attempts = " << r.num_resolve_attempts << std::endl;
	out << "crossed_path_indices = [";
	for (size_t i = 0; i < r.crossed_path_indices.size(); ++i)
		out << (i ? ", " : "") << r.crossed_path_indices[i];
	out << "]" << std::endl;
	return out;
}

/**
\brief A concise, human-readable summary of a zero-dimensional solve: how every path ended up, and
-- crucially -- whether any path was lost.

A raw solution count can lie: if a near-singular path fails to track, a genuine root is silently
missing and the count comes back short.  This report buckets every endpoint into a finite solution,
a clean divergence, or a genuine FAILURE (the tracker gave up), records the failures by their named
SuccessCode, and exposes all_paths_resolved as the one-line health check.

\see ZeroDim::Report, SummarizeSolve
*/
struct SolveReport
{
	unsigned long long num_paths_tracked = 0;    ///< total paths tracked (the Bezout / start count)
	unsigned long long num_finite_solutions = 0; ///< DISTINCT finite solutions (round of sum 1/multiplicity)
	unsigned long long num_finite_endpoints = 0; ///< raw count of finite, successful endpoints
	unsigned long long num_diverged = 0;         ///< paths ending at infinity (Success & !finite, GoingToInfinity, or SecurityMaxNormReached truncation)
	unsigned long long num_failed = 0;           ///< paths the tracker could not resolve (no solution, no clean divergence)
	unsigned long long num_singular = 0;         ///< finite solutions flagged singular (multiple / ill-conditioned)
	unsigned long long num_real = 0;             ///< finite solutions flagged real
	unsigned long long num_nonsolutions = 0;     ///< finite endpoints that are NOT solutions of the target (nonsolutions of an over-determined, squared-up system)
	std::map<SuccessCode, unsigned long long> failures_by_reason; ///< histogram of the failed paths' SuccessCodes
	double max_condition_number = 0;             ///< largest condition number among finite solutions
	unsigned max_precision_used = 0;             ///< highest working precision any path needed (digits)
	MidpathCheckReport midpath;                  ///< the path-crossing check outcome
	bool all_paths_resolved = false;             ///< num_failed==0 AND no unresolved crossings: the solve is trustworthy
};

/**
\brief Build a SolveReport from the per-endpoint final metadata and the midpath-crossing report.
Pure aggregation -- one pass, no solve state -- so it is unit-testable on synthetic metadata.
*/
template<typename ComplexT>
SolveReport SummarizeSolve(std::vector<SolutionMetaData<ComplexT>> const& metadata,
                           MidpathCheckReport const& midpath)
{
	SolveReport r;
	r.num_paths_tracked = metadata.size();
	r.midpath = midpath;
	double finite_distinct = 0;
	for (auto const& m : metadata)
	{
		if (m.max_precision_used > r.max_precision_used)
			r.max_precision_used = m.max_precision_used;

		// A path the endgame's security check truncated near infinity (SecurityMaxNormReached) is
		// diverging, not failing -- it is the security/max-norm heuristic catching a path on its way
		// to infinity, exactly like a clean GoingToInfinity.  (The endgame tests classify the two
		// together too.)  With infinite-path truncation on by default, these are expected, so they
		// must not inflate num_failed / clear all_paths_resolved.
		bool diverged = (m.endgame_success_code == SuccessCode::GoingToInfinity
		              || m.endgame_success_code == SuccessCode::SecurityMaxNormReached);
		if (m.endgame_success_code == SuccessCode::Success)
		{
			if (m.is_nonsolution)
				++r.num_nonsolutions;       // a nonsolution: finite, but not a solution of the target
			else if (m.is_finite)
			{
				++r.num_finite_endpoints;
				finite_distinct += 1.0 / m.multiplicity;
				if (m.is_singular) ++r.num_singular;
				if (m.is_real)     ++r.num_real;
				if (static_cast<double>(m.condition_number) > r.max_condition_number)
					r.max_condition_number = static_cast<double>(m.condition_number);
			}
			else
				diverged = true;
		}

		if (diverged)
			++r.num_diverged;
		else if (m.endgame_success_code != SuccessCode::Success)   // neither a solution nor a clean divergence
		{
			++r.num_failed;
			++r.failures_by_reason[m.endgame_success_code];
		}
	}
	r.num_finite_solutions = static_cast<unsigned long long>(finite_distinct + 0.5);
	r.all_paths_resolved = (r.num_failed == 0) && midpath.passed;
	return r;
}

/// \brief Stream insertion: write a human-readable summary of the solve report.
inline
std::ostream& operator<<(std::ostream & out, const SolveReport & r)
{
	out << "zero-dim solve -- " << r.num_paths_tracked << " paths tracked\n";
	out << "  finite solutions    " << r.num_finite_solutions << "   (distinct)\n";
	out << "  diverged            " << r.num_diverged << "\n";
	out << "  FAILED              " << r.num_failed;
	if (r.num_failed)
	{
		out << "    ";
		bool first = true;
		for (auto const& kv : r.failures_by_reason)
		{
			out << (first ? "" : ", ") << kv.second << " x " << kv.first;
			first = false;
		}
	}
	out << "\n  ----\n";
	out << "  singular solutions  " << r.num_singular << "\n";
	out << "  real solutions      " << r.num_real << "\n";
	if (r.num_nonsolutions)
		out << "  nonsolutions        " << r.num_nonsolutions << "   (finite, but not solutions of the target)\n";
	out << "  path crossings      " << r.midpath.num_crossings_detected
	    << (r.midpath.passed ? " (resolved)" : " (UNRESOLVED)") << "\n";
	out << "  max condition num   " << r.max_condition_number << "\n";
	out << "  max precision used  " << r.max_precision_used << " digits\n";
	out << "  all paths resolved? " << (r.all_paths_resolved ? "yes" : "NO") << "\n";
	return out;
}

/**
\brief The continuation engine -- track a set of start points through a (caller-owned) homotopy,
run the endgame, classify the endpoints, report.  See the forward-declare doc above.
*/
		template<typename TrackerType, typename EndgameType, typename SystemType>
		struct HomotopySolver :
							public virtual AnyZeroDim,
							public Observable,
							public detail::Configured<
								typename AlgoTraits< HomotopySolver<TrackerType, EndgameType, SystemType>>::NeededConfigs>
		{
			// these usings are for getters in python
			using TrackerT          = TrackerType;   ///< The path-tracker type.
			using EndgameT          = EndgameType;   ///< The endgame type.
			using SystemT           = SystemType;    ///< The system type.
			// the engine sees the start system only through the polymorphic base
			using StartSystemT       = bertini::start_system::StartSystem;  ///< The start-system type (seen polymorphically).
			using StartSystemBaseT   = bertini::start_system::StartSystem;  ///< The start-system base type.

			// This engine emits its lifecycle events on the AnyZeroDim base, so a
			// single observer type can watch any templated solver.  Accept observers
			// declared for AnyZeroDim (in addition to the exact concrete type).
			bool ObservableIsA(std::type_index t) const override
			{
				return Observable::ObservableIsA(t) || t == std::type_index(typeid(AnyZeroDim));
			}




/// a bunch of using statements to reduce typing.
			using BaseComplexT 	= typename tracking::TrackerTraits<TrackerType>::BaseComplexT;  ///< The complex number type of the tracker.
			using BaseRealT    	= typename tracking::TrackerTraits<TrackerType>::BaseRealT;     ///< The real number type of the tracker.

			using PrecisionConfig 	= typename tracking::TrackerTraits<TrackerType>::PrecisionConfig;  ///< The tracker's precision-configuration type.

			using SolnIndT 			= typename SolnCont<BaseComplexT>::size_type;  ///< Index type for solutions/paths.


			/// The Configured base storing this solver's configuration structs.
			using Config = detail::Configured<
								typename AlgoTraits<HomotopySolver<TrackerType, EndgameType, SystemType>>::NeededConfigs>;
			/// Retrieve a stored configuration struct by type (inherited from Configured).
			using Config::Get;


			using Tolerances = TolerancesConfig;          ///< Tolerances configuration type.
			using PostProcessing = PostProcessingConfig;  ///< Post-processing configuration type.
			using ZeroDimConf = ZeroDimConfig;            ///< Zero-dim configuration type.
			using AutoRetrack = AutoRetrackConfig;        ///< Auto-retrack configuration type.


			using EGBoundaryMetaDataT = EGBoundaryMetaData<BaseComplexT>;  ///< Endgame-boundary metadata type.
			using SolutionMetaDataT = SolutionMetaData<BaseComplexT>;      ///< Per-solution metadata type.


// a few more using statements

			using MidpathType = MidpathChecker<BaseRealT, BaseComplexT, EGBoundaryMetaData<BaseComplexT>>;  ///< The midpath (path-crossing) checker type.

			// --- system storage ---------------------------------------------------------------
			// The engine holds the homotopy, start system, and target by REFERENCE; the caller
			// (a user, or ZeroDimSolver) owns them and guarantees they outlive the engine.  This
			// references, not owned: the engine never forms its systems (ZeroDimSolver / the user does).
		protected:
			std::reference_wrapper<const SystemType>        target_system_;  ///< Reference to the caller-owned target system.
			std::reference_wrapper<const StartSystemBaseT>  start_system_;   ///< Reference to the caller-owned start system.
			std::reference_wrapper<const SystemType>        homotopy_;       ///< Reference to the caller-owned homotopy.

			// Re-seat the system references.  ZeroDimSolver calls this after an MPI broadcast
			// re-materializes rank 0's authoritative systems into its owned storage.
			/// \brief Re-seat the system references (used after an MPI broadcast re-materializes the owned systems).
			void ResetSystems(SystemType const& hom, StartSystemBaseT const& start, SystemType const& target)
			{
				homotopy_      = std::cref(hom);
				start_system_  = std::cref(start);
				target_system_ = std::cref(target);
			}

#ifdef BERTINI2_HAVE_MPI
			// Hook: install rank 0's authoritative systems on every rank before a distributed solve.
			// A user-supplied homotopy (plain HomotopySolver) is the user's responsibility across
			// ranks, so this is a no-op here; ZeroDimSolver, which owns its systems, overrides it.
			virtual void DistributeSystems(MPI_Comm /*comm*/) {}
#endif

		public:
			/// \return The target system being solved.
			const SystemType&       TargetSystem() const { return target_system_.get(); }
			/// \return The start system.
			const StartSystemBaseT& StartSystem()  const { return start_system_.get();  }
			/// \return The homotopy.
			const SystemType&       Homotopy()     const { return homotopy_.get();      }



/// constructors

			/**
			Construct the continuation engine over a caller-owned homotopy, start system, and target.

			- `target`:   the system the solutions satisfy at target time -- dehomogenize / residual /
			              classification reference it.
			- `start`:    supplies the start points (StartSystem::StartPoint), polymorphically.
			- `homotopy`: the system with a path variable, tracked from start to target time.
			None are owned; all three must outlive the engine.  (ZeroDimSolver builds and owns them.)
			Argument order matches the old user-homotopy constructor (target, start, homotopy).
			*/
			HomotopySolver(SystemType const& target, StartSystemBaseT const& start, SystemType const& homotopy)
			 : target_system_(std::cref(target)), start_system_(std::cref(start)), homotopy_(std::cref(homotopy)),
			   tracker_(homotopy), endgame_(tracker_)
			{
				DefaultSetup();
			}


			/**
			\brief Main Run() function provided for calling from the blackbox mode
			*/
			void Run() override
			{
				// TODO(Phase 2 Python): AnyZeroDim lacks GetSolutions() -- solutions are only
				// accessible via the concrete type's member arrays or WriteMainData().
				// Resolve when designing Python bindings; the right approach depends on whether
				// Python holds AnyZeroDim or the concrete ZeroDim<...> type directly.
#ifdef BERTINI2_HAVE_MPI
				if (parallel::Size() > 1)
				{
					RunParallel(parallel::WorldComm());
					return;
				}
#endif
				Solve();
			}

#ifdef BERTINI2_HAVE_MPI
			/**
			\brief Manager-worker parallel solve using plain C MPI.

			All ranks must have already constructed an identical ZeroDim from the same input.
			Rank 0 acts as manager; ranks 1..P-1 are workers.

			Phase 1: workers track paths to the endgame boundary.
			         Manager gathers results, runs EGBoundaryAction (midpath check).
			Phase 2: manager distributes successful paths (with boundary data) to workers.
			         Workers run the endgame, send final solutions back.
			Manager runs PostEGAction (multiplicities, same-point detection).
			*/
			void RunParallel(MPI_Comm comm)
			{
				using Result = parallel::FullPathResult<BaseComplexT>;
				using Task   = parallel::StartPointTask<BaseComplexT>;

				mpfr_free_cache(); // reproducibility: clear this rank's mpfr constant cache (see Solve())

				solutions_user_coords_fresh_ = false;

				// Each rank built its homotopy independently in the constructor -- with its OWN
				// random patch, start-system coefficients, and gamma -- so a given path index would
				// mean a different path on each worker, and the manager would gather a scrambled,
				// mostly-wrong solution set.  Make rank 0 the single authoritative source: broadcast
				// its actual target system, start system, and homotopy to every rank.  These are all
				// exact (rational coefficients + integer/rational patch), so serialization carries
				// them bit-for-bit -- there is no per-rank re-derivation that could drift, and the
				// distributed homotopy is identical to rank 0's serial one.  We also broadcast rank
				// 0's RNG seed so any randomness drawn later (e.g. in the endgame) is consistent; the
				// per-path tracking RNG is already deterministic in the path index (ReseedThisThread).
				{
					unsigned long seed = parallel::IsManager() ? GetGlobalSeed() : 0ul;
					MPI_Bcast(&seed, 1, MPI_UNSIGNED_LONG, 0, comm);
					SetGlobalSeed(seed);

					// Install rank 0's authoritative systems on every rank.  The plain engine over a
					// user-supplied homotopy does nothing here (the user owns cross-rank consistency);
					// ZeroDimSolver, which owns its systems, overrides DistributeSystems to broadcast
					// rank 0's target / start / homotopy and re-seat the engine's references.
					DistributeSystems(comm);

					num_start_points_ = StartSystem().NumStartPoints();
					GetTracker().SetSystem(Homotopy());
				}

				PreSolveChecks();
				PreSolveSetup();

				const int n_threads = parallel::WorkerThreadCount();

				// One round = dispatch a set of path indices, each executed as a WHOLE path
				// (pre-endgame + endgame) by a worker.  Round 0 is every path; later rounds re-run
				// only the crossed paths with escalated settings -- the re-track is path-independent,
				// so it rides the same demand-driven worker pool (no idle-manager bottleneck).
				if (parallel::IsManager())
				{
					// Rank 0 computes every start point once, authoritatively, and ships each to the
					// worker that tracks it.  Cached so re-track rounds reuse the identical points.
					std::vector<Vec<BaseComplexT>> start_points(num_start_points_);
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
						start_points[ii] = ComputeStartPoint(static_cast<SolnIndT>(ii));

					auto run_round = [&](std::queue<Task>& q)
					{
						parallel::RunManagerLoop<Task, Result>(comm, q,
							[this](Result const& r){ StoreFullPathResult(r); });
					};

					std::queue<Task> queue;
					for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
						queue.push(Task{ static_cast<SolnIndT>(ii), start_points[ii] });
					run_round(queue);

					// Midpath/crossing check on the collected boundary points, then bounded parallel
					// resolve rounds.  A continue/stop flag is broadcast to the workers each round so
					// they know whether to run again (and escalate their settings first).
					auto passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
					midpath_report_ = MidpathCheckReport{};
					midpath_report_.passed = passed;
					for (auto const& v : midpath_.GetCrossedPaths())
						midpath_report_.crossed_path_indices.push_back(v.index());
					midpath_report_.num_crossings_detected =
						static_cast<unsigned>(midpath_report_.crossed_path_indices.size());

					const auto max_attempts = this->template Get<ZeroDimConf>().max_num_crossed_path_resolve_attempts;
					unsigned num_resolve_attempts = 0;
					while (true)
					{
						int keep_going = (!passed && num_resolve_attempts < max_attempts) ? 1 : 0;
						MPI_Bcast(&keep_going, 1, MPI_INT, 0, comm);
						if (!keep_going)
							break;

						std::queue<Task> redo;
						for (auto const& v : midpath_.GetCrossedPaths())
							if (v.rerun())
							{
								auto idx = static_cast<SolnIndT>(v.index());
								redo.push(Task{ idx, start_points[idx] });
							}
						run_round(redo);

						passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
						++num_resolve_attempts;
					}
					midpath_report_.num_resolve_attempts = num_resolve_attempts;
					midpath_report_.passed = passed;
					if (!passed)
						std::cerr << "warning: " << midpath_report_.num_crossings_detected
						          << " path crossing(s) detected at the endgame boundary remained unresolved "
						          << "after " << num_resolve_attempts << " re-track attempt(s); "
						          << "the affected solutions may be wrong.  Consider a higher-order predictor "
						          << "or a tighter tracking tolerance.  See EndgameBoundaryMetadata()." << std::endl;

					PostEGAction();
				}
				else // worker
				{
					auto run_worker_round = [&]()
					{
						if (n_threads <= 1)
						{
							// Single-threaded worker: execute whole paths on the rank's member
							// tracker/endgame directly, from rank 0's authoritative start point.
							parallel::RunWorkerLoop<Task, Result>(comm,
								[this](Task const& task){ ExecuteOnePath(MemberDuringEGContext(), static_cast<SolnIndT>(task.path_index), task.start_point); },
								[this](Task const& task) -> Result { return PackFullPathResult(static_cast<SolnIndT>(task.path_index)); });
						}
						else
						{
							// Threaded worker: each thread owns a self-contained PathThreadState clone
							// so concurrent whole-path execution shares no mutable state.  Cloned per
							// round so it inherits any escalation (predictor bump) applied to the rank's
							// member tracker between rounds.
							auto state_factory = [this]() -> std::unique_ptr<PathThreadState>
							{
								// Clone(), not copy: System copies are SHALLOW (shared expression-tree
								// nodes whose value caches mutate on Eval); Clone() deep-copies the tree.
								// Trackers/Endgames have no default ctor, so aggregate-init via new.
								std::unique_ptr<PathThreadState> s(
									new PathThreadState{ Clone(GetTracker().GetSystem()), Clone(TargetSystem()),
									                     GetTracker(), GetEndgame(), {}, {} });
								s->tracker.SetSystem(s->sys);
								s->endgame.SetTracker(s->tracker);
								SetThreadPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
								return s;
							};

							auto track_fn = [this](std::unique_ptr<PathThreadState>& state, Task const& task) -> Result
							{
								auto idx = static_cast<SolnIndT>(task.path_index);
								ExecuteOnePath(state->Context(), idx, task.start_point);
								return PackFullPathResult(idx);
							};

							parallel::RunWorkerLoopThreaded<Task, Result>(
								comm, state_factory, track_fn, n_threads);
						}
					}; // run_worker_round

					while (true)
					{
						run_worker_round();

						int keep_going = 0;
						MPI_Bcast(&keep_going, 1, MPI_INT, 0, comm);
						if (!keep_going)
							break;

						// Match the escalation the serial RunMidpathResolution applies between
						// re-track passes: tighten midpath_retrack_tolerance_ (read by ExecuteOnePath)
						// and bump a low-order predictor on this rank's member tracker, so the next
						// round's thread clones inherit it.
						EscalateRetrackSettings();
					}
				}
			}
#endif // BERTINI2_HAVE_MPI

			// Definitions are out-of-line at the bottom of this file, after
			// output.hpp is included (avoiding a circular-include chicken-and-egg).
			void WriteMainData(std::ostream& out) const override;
			void WriteRawData(std::ostream& out)  const override;
			void WriteFiniteSolutions(std::ostream& out)      const override;
			void WriteRealFiniteSolutions(std::ostream& out)  const override;
			void WriteNonsingularSolutions(std::ostream& out) const override;
			void WriteSingularSolutions(std::ostream& out)    const override;
			void WriteRawSolutions(std::ostream& out)         const override;
			void ApplyParsedConfigs(std::string const& config_str) override;

			virtual ~HomotopySolver() = default;
/// setup functions
			// NOTE: feasibility checking (ConsistencyCheck) and start-system construction live in
			// ZeroDimSolver, which owns the systems; the engine is handed a ready homotopy.


			void DefaultSetup()
			{
				DefaultSettingsSetup();
				DefaultSystemSetup();
				DefaultTrackerSetup();
				DefaultMidpathSetup();
			}




			/// \brief Read off the number of start points the (already-built) start system will produce.
			void DefaultSystemSetup()
			{
				// The systems are already built and referenced (the engine does not form them);
				// just read off the number of start points the start system will produce.
				num_start_points_ = StartSystem().NumStartPoints(); // populate the internal variable
			}



			/**
			Fills the tolerances and retrack settings from default values.
			*/
			void DefaultSettingsSetup()
			{
				// this code can be made generic using Boost.Hana.  
				// see https://stackoverflow.com/questions/28764085/how-to-create-an-element-for-each-type-in-a-typelist-and-add-it-to-a-vector-in-c,
				// for example
				//
				// all that would need to be done is to extract a hana::tuple_t from the 
				// typelist contained in this::Config, and then hana::for_each() over it.

				this->template Set<Tolerances>(Tolerances());
				this->template Set<PostProcessing>(PostProcessing());
				// ZeroDimConfig is not templated on the complex type, so it cannot pick its own
				// per-type ambient-precision default.  Set it here, where the tracking type is known.
				ZeroDimConf zdc;
				zdc.initial_ambient_precision = DefaultInitialAmbientPrecision<BaseComplexT>();
				this->template Set<ZeroDimConf>(zdc);
				this->template Set<AutoRetrack>(AutoRetrack());
			}

			/// \brief Set the tolerance used when re-tracking paths flagged by the midpath check.
			void SetMidpathRetrackTol(NumErrorT const& rt)
			{
				midpath_retrack_tolerance_ = rt;
			}

			/// \return The tolerance used when re-tracking paths flagged by the midpath check.
			const auto& MidpathRetrackTol() const
			{
				return midpath_retrack_tolerance_;
			}

			/**
			\brief Set the precision at which start points are computed.

			Defaults to the initial ambient precision.  Set this to carry a higher precision forward
			(e.g. when one solve's output seeds the next).  Once set, PreSolveSetup will not override it.
			*/
			void SetStartPointPrecision(unsigned p)
			{
				start_point_precision_ = p;
				start_point_precision_set_by_user_ = true;
			}

			/// \return The precision at which start points are computed.
			unsigned StartPointPrecision() const
			{
				return start_point_precision_;
			}


			/**
			call this after setting up the tolerances, etc.
			*/
			void DefaultMidpathSetup()
			{
				midpath_ = MidpathType(MidPathConfig());
			}


			/// \brief Apply a configuration to the midpath (path-crossing) checker.
			void SetMidpath(MidPathConfig const& mp)
			{
				midpath_.Set(mp);
			}

			/**
			\brief Takes the default action to set up the zero dim algorithm with the default constructed tracker.

			\note should be called after the homotopy is set up, ideally.
			*/
			void DefaultTrackerSetup()
			{
				tracker_.SetSystem(Homotopy());
				tracker_.Setup(tracking::predict::DefaultPredictor(),
				              	this->template Get<Tolerances>().newton_before_endgame,
				              	this->template Get<Tolerances>().path_truncation_threshold,
								tracking::SteppingConfig(), tracking::NewtonConfig());

				tracker_.PrecisionSetup(PrecisionConfig(Homotopy()));
			}


			/**
			\brief Sets the tracker to one you supply to this function.

			Setting the tracker associates the endgame with the tracker you pass in, too.  If this is a problem, and you need this generalized so you can use a different tracker for the pre/endgame zones, please file an issue requesting this feature.

			Assumes you have done all necessary setup to it, including associating it with the homotopy for the ZeroDim algorithm.

			Again, YOU must ensure the tracker is associated with the correct homotopy
			*/
			void SetTracker(TrackerType const& new_tracker)
			{
				tracker_ = new_tracker;
				endgame_.SetTracker(tracker_);
			}

			/**
			\brief Gets the tracker

			This version gets a const reference to it.
			*/
			const TrackerType & GetTracker() const
			{
				return tracker_;
			}

			/**
			\brief Gets the tracker

			This version gets a mutable reference to it.
			*/
			TrackerType & GetTracker()
			{
				return tracker_;
			}


///  endgame specific stuff
			/**
			\brief Sets the endgame to one you supply to this function.

			Assumes you have done all necessary setup to it, including associating it with the homotopy for the ZeroDim algorithm.

			Also assumes that you have made the tracker inside the ZeroDim algorithm be the self-same tracker object as is used for the endgame you are setting here.
			*/
			void SetEndgame(EndgameType const& new_endgame)
			{
				endgame_ = new_endgame;
			}

			/**
			\brief Gets the endgame

			This version gets a `const` reference to it.
			*/
			const EndgameType & GetEndgame() const
			{
				return endgame_;
			}

			/**
			\brief Gets the endgame

			This version gets a mutable reference to it.
			*/
			EndgameType & GetEndgame()
			{
				return endgame_;
			}


			/// \brief Get a configuration struct of type \p T from the endgame.
			template<typename T>
			const T & GetFromEndgame() const
			{
				return endgame_.template Get<T>();
			}

			/// \brief Set a configuration struct of type \p T on the endgame.
			template<typename T>
			void SetToEndgame(T const& t)
			{
				endgame_.Set(t);
			}

/// the main functions



			/**
			\brief Perform the basic Zero Dim solve algorithm.

			This function iterates over all start points to the start system, tracking from each to the endgame boundary.  At the endgame boundary, path crossings are checked for.  Multiple paths which jump onto each other will not be detected, unless you are using a certified tracker, in which case this is prevented in the first place.

			Paths which have crossed are re-run, up to a certain number of times.

			Then, the points at the endgame boundary are tracked using the prescribed endgame toward the final time.

			Finally, results are post-processed.

			It is up to you to put the output somewhere.
			*/
			void Solve()
			{
				// Reproducibility: clear this thread's mpfr constant cache (pi, etc.) at the solve
				// boundary.  mpfr caches transcendental constants at the precision last requested; the
				// Cauchy endgame needs pi (roots of unity), so a high-precision solve leaves pi cached
				// at high precision.  A subsequent lower-precision solve would otherwise reuse that
				// cached value rounded down -- differing by ~1 ULP from a freshly-computed
				// low-precision pi -- making a solve's behaviour depend on what ran before it in the
				// process (non-reproducible).  Clearing here makes each solve independent of history.
				// Cheap (one constant recompute per solve); thread-local (each worker thread clears its
				// own when it starts a solve).
				mpfr_free_cache();

				solutions_user_coords_fresh_ = false;

				PreSolveChecks();

				PreSolveSetup();

				this->NotifyObservers(AlgorithmStarted<AnyZeroDim>(*this));

				// Speculative-full-path model: carry every path all the way through (pre-endgame +
				// endgame) as one unit, then detect crossings from the collected boundary points and
				// re-run any crossed path in full.  This is the same shape the distributed solve uses
				// (one worker per whole path), so serial and distributed share the per-path primitive
				// -- including computing the start point via the same ComputeStartPoint.

				// Compute every start point once, up front (deterministic per index, as the MPI
				// manager does).  The threaded dispatch ships each worker its start point, and the
				// serial path / crossed-path re-tracks reuse the identical points from here.
				std::vector<Vec<BaseComplexT>> start_points(num_start_points_);
				std::vector<SolnIndT>          all_indices(num_start_points_);
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto idx = static_cast<SolnIndT>(ii);
					start_points[ii] = ComputeStartPoint(idx);
					all_indices[ii]  = idx;
				}

				// num_threads: 0 = auto (hardware_concurrency), 1 = serial, N = N threads;
				// OMP_NUM_THREADS overrides.  n_threads <= 1 takes the pool-free serial path.
				const unsigned n_threads =
					parallel::EffectiveThreadCount(this->template Get<ZeroDimConf>().num_threads);

				if (n_threads <= 1)
				{
					for (auto idx : all_indices)
						ExecuteOnePath(MemberDuringEGContext(), idx, start_points[idx]);
				}
				else
				{
					RunPathsThreaded(all_indices, start_points, n_threads);
				}

				// Crossed-path re-tracks run on the main thread against the (escalated) member
				// tracker.  They are typically rare (often none); running them in parallel too is a
				// deferred optimization.
				RunMidpathResolution([this, &start_points](SolnIndT idx){
					ExecuteOnePath(MemberDuringEGContext(), idx, start_points[idx]);
				});

				PostEGAction();

				this->NotifyObservers(AlgorithmComplete<AnyZeroDim>(*this));
			}




			/**
			\brief Get the computed solutions, in the coordinates of the user's
			original variables.

			The stored internal solutions are dehomogenized through the target
			system, once, lazily, into a cache; repeated calls return the cached
			container by const reference.  Assumes post-solve serial access, like
			the other accessors here.

			\see SolutionsInternalCoords
			*/
			const auto& SolutionsUserCoords() const
			{
				if (!solutions_user_coords_fresh_)
				{
					solutions_user_coords_.clear();
					solutions_user_coords_.reserve(solutions_post_endgame_.size());
					for (const auto& s : solutions_post_endgame_)
					{
						// A path that failed before or during the endgame leaves its endpoint
						// slot default-constructed (zero coordinates); keep index alignment with
						// SolutionMetadata (output filters on it) by emitting an empty
						// placeholder rather than trying to dehomogenize an unset point.
						if (s.size() == 0)
							solutions_user_coords_.emplace_back();
						else
							solutions_user_coords_.push_back(this->TargetSystem().DehomogenizePoint(s));
					}
					solutions_user_coords_fresh_ = true;
				}
				return solutions_user_coords_;
			}

			/**
			\brief Get the computed solutions in the solver's internal coordinates:
			homogenized, and lying on the target system's patch.

			These are the coordinates for continuing work -- start points for further
			tracking re-using the target system's patch, refinement, etc.  Take a
			user-coordinates point back to this representation with
			System::HomogenizePoint on the target system.

			\see SolutionsUserCoords
			*/
			const auto& SolutionsInternalCoords() const
			{
				return solutions_post_endgame_;
			}

			/**
			\brief Select the solution points whose parallel metadata satisfies a predicate.

			The backing for the FiniteSolutions / RealSolutions / Singular / Nonsingular convenience
			accessors.  Returns points in user coordinates by default (pass user_coords=false for the
			solver's internal coordinates), mirroring SolutionsUserCoords / SolutionsInternalCoords.
			*/
			template<typename Pred>
			SolnCont<Vec<BaseComplexT>> SolutionsWhere(Pred pred, bool user_coords = true) const
			{
				auto const& sols = user_coords ? SolutionsUserCoords() : SolutionsInternalCoords();
				auto const& md   = SolutionMetadata();
				SolnCont<Vec<BaseComplexT>> out;
				for (size_t i = 0; i < md.size() && i < sols.size(); ++i)
					if (pred(md[i]))
						out.push_back(sols[i]);
				return out;
			}

			/**
			\brief The finite solutions: successful, finite endpoints that ARE solutions of the target
			(is_finite applies the configured endpoint_finite_threshold; nonsolutions are excluded).
			Includes singular, nonsingular, and real solutions alike.
			\see RealSolutions, SingularSolutions, NonsingularSolutions, Nonsolutions
			*/
			SolnCont<Vec<BaseComplexT>> FiniteSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success_code == SuccessCode::Success && m.is_finite && !m.is_nonsolution; }, user_coords);
			}

			/// \brief The real finite solutions (is_real applies the configured tolerance).
			SolnCont<Vec<BaseComplexT>> RealSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success_code == SuccessCode::Success && m.is_finite && !m.is_nonsolution && m.is_real; }, user_coords);
			}

			/// \brief The nonsingular finite solutions (simple, well-conditioned roots).
			SolnCont<Vec<BaseComplexT>> NonsingularSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success_code == SuccessCode::Success && m.is_finite && !m.is_nonsolution && !m.is_singular; }, user_coords);
			}

			/// \brief The singular finite solutions (multiple or ill-conditioned roots).
			SolnCont<Vec<BaseComplexT>> SingularSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success_code == SuccessCode::Success && m.is_finite && !m.is_nonsolution && m.is_singular; }, user_coords);
			}

			/**
			\brief The NONSOLUTIONS: finite, successful endpoints that are NOT solutions of the target
			system -- the extraneous nonsolutions squaring up an over-determined system introduces.
			Empty for a system solved without randomization.  These are excluded from FiniteSolutions /
			RealSolutions / etc.; a regeneration cascade reads them to discard nonsolutions.  \see is_nonsolution
			*/
			SolnCont<Vec<BaseComplexT>> Nonsolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){
					return m.endgame_success_code == SuccessCode::Success && m.is_nonsolution; }, user_coords);
			}

			/**
			\brief The solutions at infinity: endpoints not classified finite (is_finite is false) --
			the complement of FiniteSolutions within the full solution list.  These are the paths the
			endgame resolved as diverging (its GoingToInfinity / SecurityMaxNormReached verdict, or a
			successful endpoint whose dehomogenized infinity norm exceeds endpoint_finite_threshold).

			\note A path that FAILED before the endgame also leaves is_finite at its default (false) and
			so appears here; its stored point is not a meaningful solution at infinity.  Consult the
			endgame_success_code in SolutionMetadata (or Report()) to distinguish a true divergence from a
			tracking failure.  \see FiniteSolutions
			*/
			SolnCont<Vec<BaseComplexT>> InfiniteSolutions(bool user_coords = true) const
			{
				return SolutionsWhere([](auto const& m){ return !m.is_finite; }, user_coords);
			}

			/**
			\brief Get the metadat associated with the final computed solutions
			*/
			const auto& SolutionMetadata() const
			{
				return solution_final_metadata_;
			}

			/**
			\brief Get the solutions as computed at the endgame boundary
			*/
			const auto& EndgameBoundarySolutions() const
			{
				return solutions_at_endgame_boundary_;
			}

			/**
			\brief Get the report from the midpath (path-crossing) check performed at the endgame
			boundary: how many crossings were detected, which paths, how many re-track attempts were
			made, and whether the check ultimately passed.

			\see MidpathCheckReport
			*/
			const MidpathCheckReport& EndgameBoundaryMetadata() const
			{
				return midpath_report_;
			}

			/**
			\brief A concise summary of the solve: how every path ended, by category, and whether any
			path was lost.  Computed on demand from the final solution metadata and the midpath report.

			\see SolveReport
			*/
			SolveReport Report() const
			{
				return SummarizeSolve(SolutionMetadata(), EndgameBoundaryMetadata());
			}

		private:

			/**
			\brief Check that the solver functor is ready to go.
			*/
			void PreSolveChecks() const
			{
				if (num_start_points_ > solutions_at_endgame_boundary_.max_size())
					throw std::runtime_error("start system has more solutions than container for results.  I refuse to continue until this has been addressed.");
			}


			void PreSolveSetup()
			{
				// Fixed-multiple precision: FixedPrecisionConfig.precision (on the tracker) is the
				// authoritative precision for the whole solve.  Lift the tracker, the ambient/thread
				// precision, the start-point precision (via initial_ambient_precision), and the systems
				// all to that one value -- TrackerLoopInitialization requires them to agree.  (Double is
				// fixed at 16 and adaptive manages its own precision, so neither enters here.)
				if constexpr (!tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					auto const& fp = GetTracker().template Get<PrecisionConfig>();
					// Double validates (its PrecisionSetup throws on any precision other than 16);
					// fixed-multiple adopts the value into precision_.
					GetTracker().PrecisionSetup(fp);
					if constexpr (!std::is_same<BaseComplexT, complex_dbl>::value)
						if (fp.precision != 0)
						{
							ZeroDimConf zdc = this->template Get<ZeroDimConf>();
							zdc.initial_ambient_precision = fp.precision;  // -> thread + start-point precision
							this->template Set<ZeroDimConf>(zdc);
							TargetSystem().precision(fp.precision);
							Homotopy().precision(fp.precision);
						}
				}

				auto num_as_size_t = static_cast<SolnIndT>(num_start_points_);

				solution_final_metadata_.resize(num_as_size_t);
				solutions_at_endgame_boundary_.resize(num_as_size_t);
				solutions_post_endgame_.resize(num_as_size_t);

				SetMidpathRetrackTol(this->template Get<Tolerances>().newton_before_endgame);

				// Default the start-point precision to the initial ambient precision.  A caller can
				// override it via SetStartPointPrecision (e.g. to carry a higher precision forward
				// from a previous solve whose output feeds this one).
				if (!start_point_precision_set_by_user_)
					start_point_precision_ = this->template Get<ZeroDimConf>().initial_ambient_precision;
			}

			/**
			Reference-bundle over the tracker/endgame/systems/observers a single path needs.
			Lets one set of per-path execution bodies (ExecuteBeforeEG / ExecuteDuringEG) serve
			both the serial flow (members) and the distributed worker (thread-owned clones), so
			there is exactly ONE per-path code path and serial/distributed cannot diverge.
			*/
			struct BeforeEGContext
			{
				TrackerType& tracker;
				tracking::FirstPrecisionRecorder<TrackerType>&  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType>& min_max_prec;
			};
			struct DuringEGContext
			{
				// const ref: the residual/dehomogenize/precision calls are const-callable, and the
				// RefToGiven policy hands back a const target system.  A mutable thread-owned clone
				// binds here too.
				SystemType const& target_sys;
				TrackerType& tracker;
				EndgameType& endgame;
				tracking::FirstPrecisionRecorder<TrackerType>&  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType>& min_max_prec;
			};

			BeforeEGContext MemberBeforeEGContext()
			{
				return BeforeEGContext{ GetTracker(), first_prec_rec_, min_max_prec_ };
			}
			DuringEGContext MemberDuringEGContext()
			{
				return DuringEGContext{ TargetSystem(), GetTracker(), GetEndgame(), first_prec_rec_, min_max_prec_ };
			}

			/**
			\brief Compute the start point for one path, authoritatively.

			The single source of a start point.  The serial flow calls this directly; in a distributed
			solve rank 0 calls it and ships the result to the worker (the worker never derives its own),
			so start points are identical across serial and distributed runs.  Computed at
			start_point_precision_ (defaults to the initial ambient precision; settable so a chained
			solve can carry precision forward).  StartPoint() mutates the shared start system's
			expression-tree value caches, so generation is serialized.
			*/
			Vec<BaseComplexT> ComputeStartPoint(SolnIndT soln_ind)
			{
				SetThreadPrecision(start_point_precision_);
				static std::mutex start_point_mutex;
				std::lock_guard<std::mutex> lock(start_point_mutex);
				return StartSystem().template StartPoint<BaseComplexT>(soln_ind);
			}

			/**
			\brief Track one path from the start time to the endgame boundary, against `ctx`.

			The single before-endgame body, shared by the serial flow and the distributed worker.
			Uses SetThreadPrecision (thread-local) throughout, so it is correct on the main thread
			and concurrently on worker threads alike.  Writes the boundary point + the pre-endgame
			portion of the metadata.
			*/
			void ExecuteBeforeEG(BeforeEGContext ctx, SolnIndT soln_ind, Vec<BaseComplexT> const& start_point)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind));

				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					ctx.tracker.AddObserver(ctx.first_prec_rec);
					ctx.tracker.AddObserver(ctx.min_max_prec);
				}

				auto& smd = solution_final_metadata_[soln_ind];
				smd.path_index    = soln_ind;
				smd.solution_index = soln_ind;

				auto initial_prec = this->template Get<ZeroDimConf>().initial_ambient_precision;
				// SetThreadPrecision: writes thread-local only, safe from concurrent threads.
				SetThreadPrecision(initial_prec);

				// Draw the condition-number probe direction ONCE here, for the whole track of this
				// point (pre-endgame AND endgame).  We have just reseeded this thread's RNG to a value
				// determined solely by the path index, so the probe is deterministic per path and
				// identical whether the path is tracked in a serial loop or on a distributed worker --
				// regardless of how the paths were split across ranks.  It is NOT refreshed for the
				// endgame (which would churn the RNG across its sample-circle sub-tracks), so the same
				// direction persists start to finish.
				ctx.tracker.RefreshConditionDirection();

				// The times are stored precision-free (mpq_rational); materialize them as the tracking
				// complex type at the current working precision.
				BaseComplexT t_start           ( BaseRealT(this->template Get<ZeroDimConf>().start_time) );
				BaseComplexT t_endgame_boundary( BaseRealT(this->template Get<ZeroDimConf>().endgame_boundary) );

				// The start point is supplied by the caller (computed once, authoritatively, via
				// ComputeStartPoint) rather than regenerated here.  In a distributed solve rank 0
				// computes it and sends it with the task, so a worker never derives its own -- that,
				// plus authoritative pi (issue #156), is what makes the start point identical across
				// serial and distributed runs.

				// Reset to the configured initial step size for this fresh path.  In the
				// speculative-full-path model a path's endgame runs immediately before the next
				// path's pre-endgame tracking on the same tracker, and the endgame leaves
				// ReinitializeInitialStepSize(false) with a tiny step; without restoring it here the
				// next path would crawl from the start time with that tiny step (effectively a hang).
				ctx.tracker.ReinitializeInitialStepSize(true);

				// Begin every AMP path in hardware double precision; the AMP step criteria escalate to
				// multiprecision only where a path actually needs it (and the endgame, which clears
				// this override below, manages its own precision).  Total-degree start points are born
				// at LowestMultiplePrecision (the generator's arithmetic widens past the requested
				// digits, see issue #308), and initial_ambient_precision is itself a multiprecision
				// value, so without forcing double here the tracker starts -- and, in practice, stays
				// -- in multiprecision on every path, even trivially well-conditioned ones.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					ctx.tracker.SetStartPrecision(DoublePrecision());

				Vec<BaseComplexT> result;
				auto tracking_success = ctx.tracker.TrackPath(result, t_start, t_endgame_boundary, start_point);

				solutions_at_endgame_boundary_[soln_ind] =
					EGBoundaryMetaDataT({ result, tracking_success, ctx.tracker.CurrentStepsize(), ctx.tracker.CurrentPrecision() });

				// Clear the start-precision override so it does not leak into the endgame, which
				// shares this tracker instance for its sample circles.
				if constexpr (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					ctx.tracker.SetStartPrecision(std::nullopt);

				smd.pre_endgame_success_code = tracking_success;

				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (ctx.first_prec_rec.DidPrecisionIncrease())
					{
						smd.precision_changed = true;
						smd.time_of_first_prec_increase = ctx.first_prec_rec.TimeOfIncrease();
					}
					else
					ctx.tracker.RemoveObserver(ctx.first_prec_rec);
					ctx.tracker.RemoveObserver(ctx.min_max_prec);
					using std::max;
					smd.max_precision_used =
						max(smd.max_precision_used, ctx.min_max_prec.MaxPrecision());
				}
			}

			/**
			Self-contained per-thread state for executing one WHOLE path on a worker thread --
			used by both the standalone (MPI-less) threaded solve and the threaded MPI worker.
			Each std::thread owns one: a homotopy copy (tracked by `tracker`), a target-system copy
			(residual evaluation mutates System precision state), a Tracker and Endgame, and its own
			precision observers so observer-derived metadata matches serial runs.  Build a
			DuringEGContext from it to drive ExecuteOnePath.
			*/
			struct PathThreadState
			{
				System          sys;         // homotopy, tracked by `tracker`
				System          target_sys;  // for function residuals / dehomogenization
				TrackerType     tracker;
				EndgameType     endgame;
				tracking::FirstPrecisionRecorder<TrackerType>  first_prec_rec;
				tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec;

				DuringEGContext Context()
				{
					return DuringEGContext{ target_sys, tracker, endgame, first_prec_rec, min_max_prec };
				}
			};

			/**
			\brief Build the per-thread state factory: a callable () -> unique_ptr<PathThreadState>
			that each worker thread invokes once at startup.

			Clone() (deep copy), not a plain copy: System copies are SHALLOW (they share
			expression-tree nodes whose value caches mutate on Eval), so two threads sharing a System
			copy would race on those caches.  Clone() deep-copies the tree, giving each thread a fully
			independent homotopy + target system.  The tracker/endgame are copied and re-pointed at the
			thread's own systems, and the thread's mpfr precision is set.  Shared by the standalone
			threaded solve and the threaded MPI worker.
			*/
			auto MakeThreadStateFactory()
			{
				return [this]() -> std::unique_ptr<PathThreadState>
				{
					// Trackers/Endgames have no default ctor, so aggregate-init via new; {},{} default
					// the two precision recorders.
					std::unique_ptr<PathThreadState> s(
						new PathThreadState{ Clone(GetTracker().GetSystem()), Clone(TargetSystem()),
						                     GetTracker(), GetEndgame(), {}, {} });
					s->tracker.SetSystem(s->sys);       // re-point the copy at its own System
					s->endgame.SetTracker(s->tracker);  // ... and the endgame at that tracker
					SetThreadPrecision(this->template Get<ZeroDimConf>().initial_ambient_precision);
					return s;
				};
			}

			/**
			\brief Build the per-path track function: (PathThreadState&, StartPointTask) -> FullPathResult.

			Executes one whole path against the thread-local state, then packs the result.  The pack
			(read shared [idx] slots into a value) + main-thread StoreFullPathResult (install) split is
			what keeps the shared result arrays race-free across worker threads.  Shared by the
			standalone threaded solve and the threaded MPI worker.
			*/
			auto MakeThreadTrackFn()
			{
				return [this](std::unique_ptr<PathThreadState>& state,
				              parallel::StartPointTask<BaseComplexT> const& task)
				           -> parallel::FullPathResult<BaseComplexT>
				{
					auto idx = static_cast<SolnIndT>(task.path_index);
					ExecuteOnePath(state->Context(), idx, task.start_point);
					return PackFullPathResult(idx);
				};
			}

			/**
			\brief Track a set of whole paths across a shared-memory thread pool (no MPI).

			Each worker thread owns a PathThreadState clone; this (the main) thread submits all the
			given path indices, then collects exactly that many results and installs each serially via
			StoreFullPathResult.  `start_points` is indexed by global path index, so a subset of
			`indices` (e.g. crossed paths) reuses the identical authoritative start points.
			*/
			void RunPathsThreaded(std::vector<SolnIndT> const& indices,
			                      std::vector<Vec<BaseComplexT>> const& start_points,
			                      unsigned n_threads)
			{
				using Task   = parallel::StartPointTask<BaseComplexT>;
				using Result = parallel::FullPathResult<BaseComplexT>;

				auto state_factory = MakeThreadStateFactory();
				auto track_fn      = MakeThreadTrackFn();

				parallel::WorkerThreadPool<Task, Result, decltype(state_factory), decltype(track_fn)>
					pool(static_cast<int>(n_threads), state_factory, track_fn);

				for (auto idx : indices)
					pool.submit(Task{ static_cast<std::size_t>(idx), start_points[idx] });

				for (std::size_t k = 0; k < indices.size(); ++k)
					StoreFullPathResult(pool.collect());

				pool.shutdown();
			}

			/**
			\brief Apply one step of escalation to the tracker before re-tracking crossed paths.

			Two remedies are applied: (1) tighten the tracking tolerance, and (2) raise the ODE
			predictor to the default (RKF45) if a low-order predictor is in use -- a too-low-order
			predictor (notably Euler) is the most common cause of spurious boundary crossings, and
			no amount of tolerance-tightening on a first-order predictor is as effective as moving to
			a higher-order one.  The predictor bump is idempotent (once at RKF45's order it is a
			no-op).

			\note Future work: this single step wants to become a *configurable, ordered remedy list*
			(tighten tolerance, more Newton iterations, bump predictor, shrink min step size, ...),
			applied in sequence until exhausted.  That strategy abstraction is intentionally deferred;
			for now the escalation is a fixed, minimal two-remedy step.  Termination is guaranteed
			regardless, because EGBoundaryAction bounds the number of attempts.
			*/
			void EscalateRetrackSettings()
			{
				midpath_retrack_tolerance_ *= this->template Get<AutoRetrack>().midpath_decrease_tolerance_factor;
				GetTracker().SetTrackingTolerance(midpath_retrack_tolerance_);

				const auto default_predictor = tracking::predict::DefaultPredictor();
				if (tracking::predict::Order(GetTracker().GetPredictor())
				    < tracking::predict::Order(default_predictor))
					GetTracker().SetPredictor(default_predictor);
			}



			/**
			\brief Execute one whole path: start -> endgame boundary -> target, against `ctx`.

			The unit of work in the speculative-full-path model: a single worker/thread carries a
			path through both the pre-endgame tracking and the endgame, so the boundary state never
			leaves the executing context (no serialize-the-handoff, hence no precision-transfer gap).
			The pre-endgame tracking tolerance is `midpath_retrack_tolerance_` (== newton_before_endgame
			on the first pass, tightened on re-track passes by EscalateRetrackSettings); the endgame
			uses newton_during_endgame.
			*/
			void ExecuteOnePath(DuringEGContext ctx, SolnIndT soln_ind, Vec<BaseComplexT> const& start_point)
			{
				// Carry the executing tracker on the path events: &ctx.tracker is the member tracker
				// in a serial solve and the thread-local clone in a threaded one, so a meta-observer
				// attaches its per-path sub-observer to the tracker that actually runs this path.
				Observable const* exec_tracker = &ctx.tracker;

				// wall-clock the whole path (pre-endgame + endgame).  steady_clock: monotonic, safe on
				// the worker thread; written to this path's own metadata slot, so no cross-path race.
				auto const path_start_clock = std::chrono::steady_clock::now();
				auto stamp_path_time = [&]{
					solution_final_metadata_[soln_ind].path_time_seconds =
						std::chrono::duration<double>(std::chrono::steady_clock::now() - path_start_clock).count();
				};

				this->NotifyObservers(PathStarted<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind), exec_tracker));

				ctx.tracker.SetTrackingTolerance(midpath_retrack_tolerance_);
				ExecuteBeforeEG(BeforeEGContext{ ctx.tracker, ctx.first_prec_rec, ctx.min_max_prec }, soln_ind, start_point);

				if (solution_final_metadata_[soln_ind].pre_endgame_success_code != SuccessCode::Success)
				{
					stamp_path_time();
					this->NotifyObservers(PathComplete<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind), exec_tracker));
					return;
				}

				ctx.tracker.SetTrackingTolerance(this->template Get<Tolerances>().newton_during_endgame);
				ExecuteDuringEG(ctx, soln_ind);

				stamp_path_time();
				this->NotifyObservers(PathComplete<AnyZeroDim>(*this, static_cast<std::size_t>(soln_ind), exec_tracker));
			}


			/**
			\brief Detect path crossings from the collected boundary points and, for each crossed
			path, re-run it via `redo_one` (a full-path redo) with escalated settings.

			Shared by the serial flow and the distributed manager.  Populates midpath_report_.
			Bounded by max_num_crossed_path_resolve_attempts (0 = detect-and-report only), so it always
			terminates.  Unresolved crossings stay tracked but observable (report.passed == false +
			warning).
			*/
			template<typename RedoOne>
			void RunMidpathResolution(RedoOne&& redo_one)
			{
				auto passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());

				midpath_report_ = MidpathCheckReport{};
				midpath_report_.passed = passed;
				for (auto const& v : midpath_.GetCrossedPaths())
					midpath_report_.crossed_path_indices.push_back(v.index());
				midpath_report_.num_crossings_detected =
					static_cast<unsigned>(midpath_report_.crossed_path_indices.size());

				const auto max_attempts = this->template Get<ZeroDimConf>().max_num_crossed_path_resolve_attempts;
				unsigned num_resolve_attempts = 0;
				while (!passed && num_resolve_attempts < max_attempts)
				{
					EscalateRetrackSettings();
					for (auto const& v : midpath_.GetCrossedPaths())
						if (v.rerun())
							redo_one(static_cast<SolnIndT>(v.index()));
					passed = midpath_.Check(solutions_at_endgame_boundary_, StartSystem());
					++num_resolve_attempts;
				}

				midpath_report_.num_resolve_attempts = num_resolve_attempts;
				midpath_report_.passed = passed;

				if (!passed)
					std::cerr << "warning: " << midpath_report_.num_crossings_detected
					          << " path crossing(s) detected at the endgame boundary remained unresolved "
					          << "after " << num_resolve_attempts << " re-track attempt(s); "
					          << "the affected solutions may be wrong.  Consider a higher-order predictor "
					          << "or a tighter tracking tolerance.  See EndgameBoundaryMetadata()." << std::endl;
			}


			/**
			\brief Run the endgame on one path (boundary→target), against `ctx`.

			The single during-endgame body, shared by the serial flow and the distributed worker.
			Resumes at the precision the path was using at the boundary
			(`solutions_at_endgame_boundary_[idx].precision`) and writes the final solution + the
			endgame portion of the metadata.
			*/
			void ExecuteDuringEG(DuringEGContext ctx, SolnIndT soln_ind)
			{
				ReseedThisThread(static_cast<uint64_t>(soln_ind) + static_cast<uint64_t>(num_start_points_));

				auto& smd = solution_final_metadata_[soln_ind];
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (!smd.precision_changed)
						ctx.tracker.AddObserver(ctx.first_prec_rec);
					ctx.tracker.AddObserver(ctx.min_max_prec);
				}

				const auto& bdry_point = solutions_at_endgame_boundary_[soln_ind].path_point;

				ctx.tracker.SetStepSize(solutions_at_endgame_boundary_[soln_ind].last_used_stepsize);
				ctx.tracker.ReinitializeInitialStepSize(false);

				// Resume the endgame at the precision the path was actually using at the boundary,
				// rather than inferring it from the (always-multiprecision) point's mantissa, which
				// is unreliable for paths that tracked in double.
				auto start_prec = solutions_at_endgame_boundary_[soln_ind].precision;
				SetThreadPrecision(start_prec);

				ctx.endgame.SetBoundaryTime(BaseComplexT(BaseRealT(this->template Get<ZeroDimConf>().endgame_boundary)));
				ctx.endgame.SetTargetTime  (BaseComplexT(BaseRealT(this->template Get<ZeroDimConf>().target_time)));

				auto eg_success = ctx.endgame.Run(bdry_point);

				solutions_post_endgame_[soln_ind] = ctx.endgame.template FinalApproximation<BaseComplexT>();

				smd.endgame_success_code = eg_success;

				// an unsuccessful endgame has no final approximation, so the final-point-dependent
				// metadata cannot be computed.
				if (eg_success != SuccessCode::Success)
				{
					if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
					{
						ctx.tracker.RemoveObserver(ctx.first_prec_rec);
						ctx.tracker.RemoveObserver(ctx.min_max_prec);
					}
					return;
				}
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					if (!smd.precision_changed)
					{
						if (ctx.first_prec_rec.DidPrecisionIncrease())
						{
							smd.precision_changed = true;
							smd.time_of_first_prec_increase = ctx.first_prec_rec.TimeOfIncrease();
						}
						ctx.tracker.RemoveObserver(ctx.first_prec_rec);
					}
					ctx.tracker.RemoveObserver(ctx.min_max_prec);
					using std::max;
					smd.max_precision_used =
						max(smd.max_precision_used, ctx.min_max_prec.MaxPrecision());
				}
				if (tracking::TrackerTraits<TrackerType>::IsAdaptivePrec)
				{
					assert(Precision(solutions_post_endgame_[soln_ind])==Precision(ctx.endgame.template FinalApproximation<BaseComplexT>()));
					SetThreadPrecision(Precision(solutions_post_endgame_[soln_ind]));
					ctx.target_sys.precision(Precision(solutions_post_endgame_[soln_ind]));
				}
				smd.function_residual = static_cast<NumErrorT>(ctx.target_sys.Eval(solutions_post_endgame_[soln_ind]).template lpNorm<Eigen::Infinity>());
				smd.final_time_used = ctx.endgame.LatestTime();
				smd.condition_number = ctx.tracker.LatestConditionNumber();
				smd.newton_residual = ctx.tracker.LatestNormOfStep();

				smd.accuracy_estimate = ctx.endgame.ApproximateError();
				smd.accuracy_estimate_user_coords =
					static_cast<NumErrorT>( (ctx.target_sys.DehomogenizePoint(solutions_post_endgame_[soln_ind]) -
					ctx.target_sys.DehomogenizePoint(ctx.endgame.template PreviousApproximation<BaseComplexT>())).template lpNorm<Eigen::Infinity>() );
				smd.cycle_num = ctx.endgame.CycleNumber();

				// Honest precision/accuracy of this solution, as digit counts.  precision_digits is the
				// precision the endgame actually finished in -- which is exactly the precision of the
				// stored solution point (DoublePrecision() for a fast-lane path, the mpfr precision for an
				// escalated one).  accuracy_digits is how many of those digits the convergence agreement
				// supports: floor(-log10(accuracy_estimate)), clamped to [0, precision_digits].
				smd.precision_digits = static_cast<unsigned>(Precision(solutions_post_endgame_[soln_ind]));
				{
					using std::log10; using std::floor; using std::min; using std::max;
					double err = static_cast<double>(smd.accuracy_estimate);
					unsigned acc = (err > 0.0)
						? static_cast<unsigned>(max(0.0, floor(-log10(err))))
						: smd.precision_digits;
					smd.accuracy_digits = min(acc, smd.precision_digits);
				}
			}


		protected:
			// virtual so ZeroDimSolver can append its extraneous-solution filter after the engine's
			// classification (it overrides this to call the base, then filter against the original
			// over-determined system).  Everything from here down -- the classification helpers, the
			// pack/store helpers, and the data members -- is protected so the derived algorithm can
			// read the per-endpoint metadata and endpoints it needs for that filter.
			/// \brief Action taken after the endgame: compute post-track metadata (ZeroDimSolver overrides to also filter extraneous solutions).
			virtual void PostEGAction()
			{
				ComputePostTrackMetadata();
			}



			/**
			\brief Populates the result_ member, with the post-processed results of the zero dim solve.

			output: result_ member variable.
			*/
			void ComputePostTrackMetadata()
			{
				ClassifyFiniteAndReal();
				ComputeMultiplicities();
				ClassifySingular();
			}

			/**
			\brief Classify each successful endpoint as finite/infinite, and (when finite) real.

			Uses the same dehomogenize-then-infinity-norm measurement the endgames use for their
			`Security::max_norm` divergence test (`System::InfinityNormOfDehomogenized`), compared
			against `endpoint_finite_threshold` -- so the metadata can never contradict the
			endgame's own verdict.  Endpoints the endgame already flagged as diverging
			(`GoingToInfinity` / `SecurityMaxNormReached`) are taken as infinite without
			recomputation.  A finite endpoint is real if the infinity norm of the imaginary parts
			of its dehomogenized coordinates is below `real_threshold`.
			*/
			void ClassifyFiniteAndReal()
			{
				using std::abs; using std::imag;
				const auto& post = this->template Get<PostProcessing>();

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto& smd = solution_final_metadata_[ii];

					if (smd.endgame_success_code==SuccessCode::GoingToInfinity ||
					    smd.endgame_success_code==SuccessCode::SecurityMaxNormReached)
					{
						smd.is_finite = false; // the endgame already decided this path diverges
						continue;
					}
					if (smd.endgame_success_code!=SuccessCode::Success)
						continue;              // failed otherwise: leave defaults (not finite/real/singular)

					auto user_pt = this->TargetSystem().DehomogenizePoint(solutions_post_endgame_[ii]);

					smd.is_finite =
						static_cast<NumErrorT>(user_pt.template lpNorm<Eigen::Infinity>()) <= post.endpoint_finite_threshold;

					if (smd.is_finite)
					{
						NumErrorT max_imag{0};
						for (Eigen::Index k{0}; k < user_pt.size(); ++k)
						{
							NumErrorT a = static_cast<NumErrorT>(abs(imag(user_pt(k))));
							if (a > max_imag) max_imag = a;
						}
						smd.is_real = max_imag < post.real_threshold;
					}
				}
			}

			/**
			\brief Cluster identical endpoints to compute multiplicities.

			Two endpoints are the same point when the infinity norm of the difference of their
			*dehomogenized* coordinates is below `final_tolerance * same_point_tolerance_multiplier`.
			Comparing dehomogenized (user) coordinates -- not the internal homogenized, on-patch
			coordinates, which carry the homogenizing variable and patch scaling -- is essential;
			the infinity norm matches the convergence norm the endgames use.  Only finite,
			successful endpoints are clustered (an at-infinity endpoint has no meaningful
			dehomogenized coordinates to compare).
			*/
			void ComputeMultiplicities()
			{
				const NumErrorT same_tol =
					this->template Get<Tolerances>().final_tolerance *
					this->template Get<PostProcessing>().same_point_tolerance_multiplier;

				std::vector<Vec<BaseComplexT>> user_pts(num_start_points_);
				std::vector<char> eligible(num_start_points_, 0);
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto const& smd = solution_final_metadata_[ii];
					if (smd.endgame_success_code==SuccessCode::Success && smd.is_finite)
					{
						user_pts[ii] = this->TargetSystem().DehomogenizePoint(solutions_post_endgame_[ii]);
						eligible[ii] = 1;
					}
				}

				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					if (!eligible[ii]) continue;
					for (decltype(num_start_points_) jj{ii+1}; jj < num_start_points_; ++jj)
					{
						if (!eligible[jj]) continue;
						if ( static_cast<NumErrorT>((user_pts[ii] - user_pts[jj]).template lpNorm<Eigen::Infinity>()) < same_tol )
						{
							++solution_final_metadata_[ii].multiplicity;
							++solution_final_metadata_[jj].multiplicity;
							// jj coincides with the earlier ii, so it is a duplicate, not the
							// representative; the lowest-index member of a cluster is never marked
							// here, so it remains the single representative.
							solution_final_metadata_[jj].multiplicity_representative = false;
						}
					}
				}
			}

			/**
			\brief Classify each successful endpoint as singular or not.

			Matching Bertini 1: an endpoint is singular if it is the endpoint of multiple paths
			(multiplicity > 1), or if the approximation of its condition number (spectral norm, as
			estimated by the tracker) exceeds `condition_number_threshold`.
			*/
			void ClassifySingular()
			{
				const NumErrorT cond_threshold = this->template Get<PostProcessing>().condition_number_threshold;
				for (decltype(num_start_points_) ii{0}; ii < num_start_points_; ++ii)
				{
					auto& smd = solution_final_metadata_[ii];
					if (smd.endgame_success_code!=SuccessCode::Success)
						continue;
					smd.is_singular = (smd.multiplicity > 1) || (smd.condition_number > cond_threshold);
				}
			}



		///////
		//	MPI pack/store helpers (only compiled when BERTINI2_HAVE_MPI is defined)
		///////

			// Pack everything computed for one whole path into a single result value.  Used as the
			// thread-safe handoff in BOTH the MPI protocol (serialized over the wire) and the
			// standalone threaded solve (a worker thread packs, the main thread installs via
			// StoreFullPathResult) -- the compute/install split is what keeps shared result arrays
			// race-free when several threads run distinct paths.
			/// \brief Pack everything computed for one whole path into a single result value (thread/MPI handoff).
			/// \param idx The path index to pack.
			/// \return The packed whole-path result.
			parallel::FullPathResult<BaseComplexT> PackFullPathResult(SolnIndT idx) const
			{
				parallel::FullPathResult<BaseComplexT> r;
				r.path_index          = idx;
				r.pre_endgame_success_code = solutions_at_endgame_boundary_[idx].success_code;
				r.boundary_point      = solutions_at_endgame_boundary_[idx].path_point;
				r.boundary_stepsize   = solutions_at_endgame_boundary_[idx].last_used_stepsize;
				r.boundary_precision  = solutions_at_endgame_boundary_[idx].precision;

				auto const& smd = solution_final_metadata_[idx];
				r.endgame_success_code   = smd.endgame_success_code;
				r.solution    = solutions_post_endgame_[idx];
				r.function_residual = smd.function_residual;
				r.condition_number  = smd.condition_number;
				r.newton_residual   = smd.newton_residual;
				r.final_time_used   = smd.final_time_used;
				r.accuracy_estimate = smd.accuracy_estimate;
				r.accuracy_estimate_user_coords = smd.accuracy_estimate_user_coords;
				r.cycle_num         = smd.cycle_num;
				r.precision_digits = smd.precision_digits;
				r.accuracy_digits  = smd.accuracy_digits;
				r.precision_changed = smd.precision_changed;
				r.time_of_first_prec_increase = smd.time_of_first_prec_increase;
				r.max_precision_used = smd.max_precision_used;
				r.path_time_seconds = smd.path_time_seconds;
				return r;
			}

			// Install a whole-path result on the manager.  Overwrites cleanly, so re-dispatching a
			// crossed path in a later resolve round simply replaces its earlier (crossed) result.
			/// \brief Install a whole-path result (from PackFullPathResult) into the shared result arrays.
			/// \param r The packed whole-path result to store.
			void StoreFullPathResult(parallel::FullPathResult<BaseComplexT> const& r)
			{
				auto idx = static_cast<SolnIndT>(r.path_index);
				solutions_at_endgame_boundary_[idx] =
					EGBoundaryMetaDataT{ r.boundary_point, r.pre_endgame_success_code, r.boundary_stepsize, r.boundary_precision };
				solutions_post_endgame_[idx] = r.solution;

				auto& smd = solution_final_metadata_[idx];
				smd.path_index          = idx;
				smd.solution_index      = idx;
				smd.pre_endgame_success_code = r.pre_endgame_success_code;
				smd.endgame_success_code     = r.endgame_success_code;
				smd.function_residual   = r.function_residual;
				smd.condition_number    = r.condition_number;
				smd.newton_residual     = r.newton_residual;
				smd.final_time_used     = r.final_time_used;
				smd.accuracy_estimate   = r.accuracy_estimate;
				smd.accuracy_estimate_user_coords = r.accuracy_estimate_user_coords;
				smd.cycle_num           = r.cycle_num;
				smd.precision_digits  = r.precision_digits;
				smd.accuracy_digits   = r.accuracy_digits;
				smd.precision_changed   = r.precision_changed;
				smd.time_of_first_prec_increase = r.time_of_first_prec_increase;
				smd.max_precision_used  = r.max_precision_used;
				smd.path_time_seconds   = r.path_time_seconds;
			}


		///////
		//	private data members
		///////

			unsigned long long num_start_points_;  ///< Number of start points the start system produces.
			NumErrorT midpath_retrack_tolerance_;  ///< Tolerance used when re-tracking paths flagged by the midpath check.
			MidpathCheckReport midpath_report_; ///< populated by EGBoundaryAction; exposed via EndgameBoundaryMetadata()
			unsigned start_point_precision_ = DoublePrecision(); ///< precision at which start points are computed; defaults to initial ambient precision (see PreSolveSetup)
			bool start_point_precision_set_by_user_ = false;  ///< Whether the start-point precision was set explicitly by the user.


			/// observers used during tracking
			// i feel like these should be factored out into some policy class which prescribes how they are used, so that the actions taken are customizable.
			tracking::FirstPrecisionRecorder<TrackerType> first_prec_rec_;   ///< Records the time of the first precision increase on each path.
			tracking::MinMaxPrecisionRecorder<TrackerType> min_max_prec_;    ///< Records the min/max precision used on each path.


			/// function objects used during the algorithm
			TrackerType tracker_;    ///< The path tracker.
			EndgameType endgame_;    ///< The endgame.
			MidpathType midpath_;    ///< The midpath (path-crossing) checker.



			/// computed data
			SolnCont< EGBoundaryMetaDataT > solutions_at_endgame_boundary_; ///< Per-path state captured at the endgame boundary.
			SolnCont<Vec<BaseComplexT> > solutions_post_endgame_;          ///< Solution points after the endgame (internal coordinates).
			mutable SolnCont<Vec<BaseComplexT> > solutions_user_coords_; ///< lazy cache for SolutionsUserCoords; serial post-solve access assumed
			mutable bool solutions_user_coords_fresh_ = false;            ///< Whether the user-coordinate cache is up to date.
			SolnCont<SolutionMetaDataT> solution_final_metadata_;        ///< Final per-solution metadata.


		}; // struct HomotopySolver


		// impl namespace deliberately NOT named `detail`: bertini::detail (TypeList, Configured)
		// is referenced unqualified throughout this file, and a bertini::algorithm::detail would
		// shadow it for any code defined after this point (e.g. ApplyParsedConfigs).
		namespace zero_dim_detail {

		/**
		\brief Owns and builds the target / start system / homotopy for a zero-dim solve.

		Constructed as the FIRST base of ZeroDimSolver -- before the HomotopySolver engine base, which
		holds references into these systems -- so the homotopy is fully formed before the engine reads
		it.  This holds the system-building half of the algorithm, plus the feasibility / rank checks the
		algorithm owns.  Members carry an `owned_` prefix so they do not collide with the engine base's
		reference members.
		*/
		template<typename SystemType>
		struct OwnedHomotopy
		{
			using StartSystemBaseT = bertini::start_system::StartSystem;                  ///< The start-system base type.
			using FactoryT = bertini::start_system::StartSystemFactory<SystemType>;       ///< The start-system factory type.

			/// \brief Build and own the target, start system, and homotopy from a target system and factory.
			OwnedHomotopy(SystemType const& target, FactoryT factory, std::string const& path_variable_name)
			 : owned_target_(Clone(target)), owned_factory_(std::move(factory))
			{
				ConsistencyCheck();              // feasibility (no path var; not under-constrained; polynomial)
				SquareUp();                      // randomize an over-determined system down to square
				RankCheck();                     // square, but can isolated solutions even exist?
				PrepareTarget(owned_target_);    // homogenize + auto-patch the (now square) target
				owned_start_    = owned_factory_(owned_target_);                         // start system over the prepared target
				// Treat the configured name as a BASE: mangle it only if it would collide with a
				// user variable, so the injected path variable can never overlap a user symbol.
				owned_homotopy_ = MakeHomotopy(owned_target_, *owned_start_,
				                               UniquePathVariableName(owned_target_, path_variable_name));
			}

			/// \return The built (cloned, prepared) target system.
			SystemType const&       BuiltTarget()   const { return owned_target_;   }
			/// \return The built start system.
			StartSystemBaseT const& BuiltStart()    const { return *owned_start_;   }
			/// \return The built homotopy.
			SystemType const&       BuiltHomotopy() const { return owned_homotopy_; }

		protected:
			SystemType                        owned_target_;    ///< The owned (cloned, prepared) target system.
			std::shared_ptr<StartSystemBaseT> owned_start_;     ///< The owned start system.
			SystemType                        owned_homotopy_;  ///< The owned homotopy.
			FactoryT                          owned_factory_;   ///< The start-system factory used to build the start system.

			// When the user's system was over-determined, SquareUp randomizes it down to square for
			// tracking and keeps the ORIGINAL (natural, un-homogenized) system here so ZeroDimSolver
			// can discard the extraneous solutions the squaring introduces.  was_randomized_ gates
			// that filter; randomization_matrix_ records the (exact) coefficient matrix used.
			bool                  was_randomized_ = false;       ///< Whether an over-determined target was squared up by randomization.
			SystemType            original_natural_target_;      ///< The original (natural, un-homogenized) target, kept to discard extraneous solutions.
			Mat<complex_mp>       randomization_matrix_;         ///< The exact coefficient matrix used to square up the system.

			/**
			\brief Reject targets that cannot have isolated solutions for structural reasons.
			*/
			void ConsistencyCheck() const
			{
				if (owned_target_.HavePathVariable())
					throw std::runtime_error("unable to perform zero dim solve on target system -- has path variable, use a HomotopySolver instead.");

				// A square zero-dim system needs one equation per dimension.  Each projective
				// (homogeneous) variable group of size k spans P^{k-1}: its k coordinates carry only
				// k-1 dimensions because scale is free, so it needs one fewer equation than it has
				// variables.  Subtract that free scale per projective group before comparing.  An
				// UNDER-determined system has a positive-dimensional solution set, so no isolated
				// solutions to compute -- raise a helpful error rather than tracking garbage.
				if (owned_target_.NumVariables() - owned_target_.NumHomVariableGroups() > owned_target_.NumTotalFunctions())
					throw std::runtime_error("unable to perform zero dim solve on target system -- it is under-determined (fewer equations than variables), so its solution set is positive-dimensional, not zero-dimensional.  ZeroDimSolver computes isolated solutions only; add equations, or use a positive-dimensional method.");

				if (!owned_target_.IsPolynomial())
					throw std::runtime_error("unable to perform zero dim solve on target system -- system is non-polynomial, use a HomotopySolver instead.");
			}

			/**
			\brief If the target is over-determined (more equations than the affine dimension), replace
			it with n generic combinations (System::Randomize) so a start system can track it.  The
			randomized system's isolated solutions contain the original's plus extraneous ones; the
			original system is kept (original_natural_target_) so ZeroDimSolver can filter those out.
			*/
			void SquareUp()
			{
				auto const dimension = owned_target_.NumVariables() - owned_target_.NumHomVariableGroups();
				if (owned_target_.NumTotalFunctions() <= dimension)
					return; // already square (or under-determined, which ConsistencyCheck already rejected)

				original_natural_target_ = Clone(owned_target_);   // keep the original N-function system
				owned_target_            = owned_target_.Randomize();
				randomization_matrix_    = owned_target_.RandomizationMatrix();
				was_randomized_          = true;
			}

			/**
			\brief Reject a square target whose solution set is still positive-dimensional.

			ConsistencyCheck only counts equations; a system can be square yet positive-dimensional
			(e.g. a repeated equation).  At a generic point a zero-dimensional system has a full-rank
			n x n Jacobian, so a rank-deficient Jacobian there means no isolated solutions.  Runs only
			when the (natural, pre-homogenize) system is square; projective/structured inputs whose
			Jacobian is not n x n are left to ConsistencyCheck.
			*/
			void RankCheck() const
			{
				if (owned_target_.HavePathVariable())
					return;
				auto const n = static_cast<Eigen::Index>(owned_target_.NumVariables());
				if (static_cast<Eigen::Index>(owned_target_.NumTotalFunctions()) != n || n == 0)
					return; // not a square map: ConsistencyCheck governs feasibility here

				// a generic complex sample point: a zero-dimensional variety misses it almost surely,
				// where the Jacobian attains its generic (full) rank.
				Vec<complex_dbl> pt = Vec<complex_dbl>::Random(n);
				Mat<complex_dbl> J  = owned_target_.template Jacobian<complex_dbl>(pt);
				Eigen::FullPivLU<Mat<complex_dbl>> lu(J);
				lu.setThreshold(1e-10);
				if (lu.rank() < n)
					throw std::runtime_error("unable to perform zero dim solve on target system -- the Jacobian is rank-deficient at a generic point, so the solution set is positive-dimensional, not zero-dimensional.  ZeroDimSolver computes isolated solutions only.");
			}

			/// \brief Homogenize and auto-patch the (already square) target system in place.
			/// \param target The target system to prepare.
			static void PrepareTarget(SystemType& target)
			{
				target.Homogenize(); // work over projective coordinates
				target.AutoPatch();  // then patch if needed
			}
		};

		} // ns zero_dim_detail


		/**
		\brief The zero-dimensional solve algorithm.

		Owns its systems: it clones the user's target, homogenizes/patches it, builds a start system
		(via the injected factory) and a homotopy, then drives the continuation engine.  It IS-A
		HomotopySolver -- the engine machinery (tracking, endgame, crossing resolution, classification,
		reporting) is reused, not duplicated -- with an OwnedHomotopy base supplying the systems the
		engine references.
		*/
		template<typename TrackerType, typename EndgameType, typename SystemType>
		struct ZeroDimSolver :
			private zero_dim_detail::OwnedHomotopy<SystemType>,
			public  HomotopySolver<TrackerType, EndgameType, SystemType>
		{
			using OwnedT   = zero_dim_detail::OwnedHomotopy<SystemType>;            ///< The system-owning base.
			using EngineT  = HomotopySolver<TrackerType, EndgameType, SystemType>;  ///< The continuation-engine base.
			using FactoryT = typename OwnedT::FactoryT;  ///< The start-system factory type.

			/**
			Build the start system + homotopy from `target`, then construct the engine over them.
			Base initialization order is declaration order: OwnedHomotopy (which builds) runs before
			the HomotopySolver engine (which references the freshly built systems).  The factory
			defaults to TotalDegreeBinomial (the interim default start system).
			*/
			ZeroDimSolver(SystemType const& target,
			              FactoryT factory = bertini::start_system::MakeStartFactory<bertini::start_system::TotalDegreeBinomial, SystemType>())
			 : OwnedT(target, std::move(factory), ZeroDimConfig{}.path_variable_name),
			   EngineT(OwnedT::BuiltTarget(), OwnedT::BuiltStart(), OwnedT::BuiltHomotopy())
			{}

			/// \brief Whether the supplied system was over-determined and squared-up by randomization.
			bool WasRandomized() const { return this->was_randomized_; }
			/// \brief The randomization matrix used to square up (empty if the system was already square).
			Mat<complex_mp> const& RandomizationMatrix() const { return this->randomization_matrix_; }

		protected:
			/**
			\brief After the engine classifies the endpoints, discard the extraneous solutions that
			squaring an over-determined system introduces.

			Squaring replaces N equations by n generic combinations, so the square system's isolated
			solutions are the genuine ones PLUS spurious points that satisfy the combinations but not
			the original system.  Re-evaluate the ORIGINAL (un-randomized) system at each finite
			endpoint; a residual above a solve-accuracy threshold flags the point `is_nonsolution`
			(it stays geometrically finite, but drops out of FiniteSolutions / Real / Singular and is
			surfaced by Nonsolutions()), and the reported function_residual is updated to that
			meaningful value.  A no-op unless the system was squared up.
			*/
			void PostEGAction() override
			{
				EngineT::PostEGAction();   // the engine's finite/real/multiplicity/singular classification

				if (!this->was_randomized_)
					return;

				// genuine roots satisfy the original system to ~solve accuracy; spurious ones miss it
				// by O(1).  A generous multiple of the endgame tolerance separates the two cleanly.
				const double threshold =
					1e3 * static_cast<double>(this->template Get<TolerancesConfig>().newton_during_endgame);

				for (decltype(this->num_start_points_) ii{0}; ii < this->num_start_points_; ++ii)
				{
					auto& smd = this->solution_final_metadata_[ii];
					if (smd.endgame_success_code != SuccessCode::Success || !smd.is_finite)
						continue;

					auto user_pt = this->TargetSystem().DehomogenizePoint(this->solutions_post_endgame_[ii]);
					// Evaluate the original system in DOUBLE precision: the filter only needs to tell a
					// genuine root (tiny residual) from an extraneous one (O(1)), and a double residual
					// does that robustly without a precision mismatch between the (possibly AMP/high
					// precision) solution point and the original system's working precision.
					Vec<complex_dbl> pt_d(user_pt.size());
					for (Eigen::Index k = 0; k < user_pt.size(); ++k)
						pt_d(k) = complex_dbl(user_pt(k));
					auto residual = static_cast<NumErrorT>(
						this->original_natural_target_.template Eval<complex_dbl>(pt_d).template lpNorm<Eigen::Infinity>());
					smd.function_residual = residual;            // report residual against the ORIGINAL system
					if (static_cast<double>(residual) > threshold)
						smd.is_nonsolution = true;               // extraneous: finite, but not a solution of the original system
				}
			}

#ifdef BERTINI2_HAVE_MPI
			// Broadcast rank 0's authoritative owned systems to every rank, then re-seat the engine's
			// references and let the caller (RunParallel) re-point the tracker at the homotopy.
			void DistributeSystems(MPI_Comm comm) override
			{
				parallel::mpi_broadcast_serialized(comm, this->owned_target_,   0);
				// the OWNING shared_ptr<StartSystem> carries the concrete derived type polymorphically
				// (BOOST_CLASS_EXPORT); a base reference would serialize only the System slice.
				parallel::mpi_broadcast_serialized(comm, this->owned_start_,    0);
				parallel::mpi_broadcast_serialized(comm, this->owned_homotopy_, 0);
				this->ResetSystems(this->owned_homotopy_, *this->owned_start_, this->owned_target_);
			}
#endif
		};

	} // ns algo

} // ns bertini

// Include output formatters after ZeroDim is fully defined.
// output.hpp includes zero_dim_solve.hpp, so #pragma once prevents re-inclusion
// and the circular dependency is resolved.  Any TU that gets zero_dim_solve.hpp
// therefore also gets output.hpp, making WriteMainData/WriteRawData instantiable.
#include "bertini2/nag_algorithms/output.hpp"

namespace bertini {
namespace algorithm {

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteMainData(std::ostream& out) const
{
	output::Classic<HomotopySolver>::MainData(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteRawData(std::ostream& out) const
{
	output::Classic<HomotopySolver>::RawData(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteFiniteSolutions(std::ostream& out) const
{
	output::Classic<HomotopySolver>::FiniteSolutions(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteRealFiniteSolutions(std::ostream& out) const
{
	output::Classic<HomotopySolver>::RealFiniteSolutions(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteNonsingularSolutions(std::ostream& out) const
{
	output::Classic<HomotopySolver>::NonsingularSolutions(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteSingularSolutions(std::ostream& out) const
{
	output::Classic<HomotopySolver>::SingularSolutions(out, *this);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
inline void
HomotopySolver<TrackerType,EndgameType,SystemType>::WriteRawSolutions(std::ostream& out) const
{
	output::Classic<HomotopySolver>::RawSolutions(out, *this);
}

} // ns algorithm
} // ns bertini

// Include config parsers after ZeroDim is fully defined, then provide the
// out-of-line definition of ApplyParsedConfigs.  Parsing headers do not include
// zero_dim_solve.hpp, so there is no circular dependency here.
#include "bertini2/io/parsing/settings_parsers.hpp"

namespace bertini {
namespace algorithm {

/// \brief Inject every element of a tuple of config structs into a Configured-derived target via Set<T>.
/// \tparam Target The Configured-derived target type.
/// \tparam Ts The config struct types in the tuple.
/// \param target The target to set the configs on.
/// \param t The tuple of config structs to inject.
// Injects every element of a std::tuple<Ts...> into `target` via Set<T>.
// Works for any Configured<>-derived target whose typelist contains all Ts.
template<typename Target, typename... Ts>
void InjectParsedTuple(Target& target, std::tuple<Ts...> const& t) {
	(target.template Set<Ts>(std::get<Ts>(t)), ...);
}

template<typename TrackerType, typename EndgameType, typename SystemType>
void
HomotopySolver<TrackerType,EndgameType,SystemType>
    ::ApplyParsedConfigs(std::string const& config_str)
{
	using namespace parsing::classic;

	// 1. solver-owned configs (Tolerances, PostProcessing, ZeroDimConf, AutoRetrack)
	using ZDConfs = typename Config::UsedConfigs;
	auto zd = ConfigParser<ZDConfs>::Parse(config_str);
	InjectParsedTuple(*this, zd);
	// The path variable name comes from ZeroDimConf, so the homotopy must be
	// re-formed if the parsed name differs from the one used at construction.
	// Only re-run setup in that case: re-running unconditionally re-randomizes
	// gamma and the start system for no reason, and (before Homogenize() was made
	// idempotent) re-corrupted the prepared target system.
	// Compare against the collision-free name actually used at construction (the configured
	// name is a base that UniquePathVariableName may have mangled), so an unchanged config
	// does not spuriously look like a rename and force a rebuild.
	auto const effective_path_variable_name =
	    UniquePathVariableName(this->TargetSystem(), this->template Get<ZeroDimConf>().path_variable_name);
	if (!Homotopy().HavePathVariable()
	    || Homotopy().GetPathVariable()->name() != effective_path_variable_name)
	{
		DefaultSystemSetup();
	}

	// 2. Tracker — uses positional Setup() rather than Set<T>, so handle explicitly.
	using TkConfs = detail::TypeList<
	    tracking::SteppingConfig,
	    tracking::NewtonConfig,
	    tracking::Predictor>;
	auto tk = ConfigParser<TkConfs>::Parse(config_str);
	tracker_.Setup(
	    std::get<tracking::Predictor>(tk),
	    this->template Get<Tolerances>().newton_before_endgame,
	    this->template Get<Tolerances>().path_truncation_threshold,
	    std::get<tracking::SteppingConfig>(tk),
	    std::get<tracking::NewtonConfig>(tk));
	tracker_.PrecisionSetup(PrecisionConfig(Homotopy()));
	endgame_.SetTracker(tracker_); // keep endgame's tracker ref consistent after reconfigure

	// 3. Endgame-owned configs — endgame_ is the concrete EndgameType deriving from
	//    Configured<its_configs>, so Set<T> is available directly.
	//    AlgoTraits<EndgameType>::NeededConfigs is public (unlike EndgameType::Configs).
	using EGConfs = typename endgame::AlgoTraits<EndgameType>::NeededConfigs;
	auto eg = ConfigParser<EGConfs>::Parse(config_str);
	InjectParsedTuple(endgame_, eg);

	// 4. Midpath config
	SetMidpath(ConfigParser<MidPathConfig>::Parse(config_str));
}

} // ns algorithm
} // ns bertini


// Explicit instantiation declarations — suppress re-instantiation of the production solver types in
// every including TU.  Six ZeroDimSolver combos (Tracker x Endgame) cover EVERY clone-owned start
// system (the start system is held polymorphically), plus six HomotopySolver combos for user
// homotopies.  Definitions in core/src/eti/zero_dim_eti.cpp + zero_dim_blackbox_eti.cpp; see ADR-0014.
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"

namespace bertini{ namespace algorithm{

extern template struct ZeroDimSolver<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::PSEG,     System>;
extern template struct ZeroDimSolver<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::Cauchy,   System>;
extern template struct ZeroDimSolver<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::PSEG,   System>;
extern template struct ZeroDimSolver<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::Cauchy, System>;
extern template struct ZeroDimSolver<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::PSEG,                 System>;
extern template struct ZeroDimSolver<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::Cauchy,               System>;

// user homotopies: the engine over caller-owned systems.  definitions in zero_dim_blackbox_eti.cpp
extern template struct HomotopySolver<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::PSEG,     System>;
extern template struct HomotopySolver<tracking::DoublePrecisionTracker,   typename endgame::EndgameSelector<tracking::DoublePrecisionTracker>::Cauchy,   System>;
extern template struct HomotopySolver<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::PSEG,   System>;
extern template struct HomotopySolver<tracking::MultiplePrecisionTracker, typename endgame::EndgameSelector<tracking::MultiplePrecisionTracker>::Cauchy, System>;
extern template struct HomotopySolver<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::PSEG,                 System>;
extern template struct HomotopySolver<tracking::AMPTracker,               typename endgame::EndgameSelector<tracking::AMPTracker>::Cauchy,               System>;

}} // namespaces
