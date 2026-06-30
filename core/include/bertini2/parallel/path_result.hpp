//This file is part of Bertini 2.
//
//bertini2/parallel/path_result.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/path_result.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/path_result.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/path_result.hpp

\brief Result struct used in the manager-worker MPI protocol for ZeroDim.

In the unified speculative-full-path model a worker executes one WHOLE path (pre-endgame +
endgame), so the protocol is a single task/result pair:
  task   = SolnIndT path index
  result = FullPathResult<ComplexT>  (boundary data + final solution + endgame metadata)

FullPathResult carries Boost.Serialization support so Boost.MPI can transmit it via the existing
serialization for Eigen vectors and arbitrary-precision types.
*/

#pragma once

#include <deque>
#include <limits>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"   // defines Vec<T> and sets up Eigen serialization plugin
#include "bertini2/mpfr_extensions.hpp"    // Boost.Serialization for real_mp, mpc_complex
#include "bertini2/common/config.hpp"      // SuccessCode (needs Vec<T> and <deque> first)

#include <boost/serialization/vector.hpp>

namespace bertini {
namespace parallel {

constexpr int TAG_WORK_ITEM = 1;  ///< MPI message tag: manager -> worker work item.
constexpr int TAG_RESULT    = 2;  ///< MPI message tag: worker -> manager result.
constexpr int TAG_CAPACITY  = 3;  ///< MPI message tag: worker -> manager, max tasks in flight on that rank.



/**
\brief Work item for the speculative-full-path model: a path index plus its start point.

Rank 0 computes start points authoritatively and ships them to workers, so a worker never derives
its own (that, with authoritative pi, is what keeps start points -- and hence whole tracks --
identical across serial and distributed runs).
*/
template<typename ComplexT>
struct StartPointTask
{
	using SolnIndT = std::size_t;  ///< The path-index type.

	SolnIndT      path_index = std::numeric_limits<SolnIndT>::max();  ///< The path index (max value marks a sentinel).
	Vec<ComplexT> start_point;  ///< The authoritative start point for this path.

	/// \brief Query whether this is the sentinel "no more work" task.
	bool is_sentinel() const
	{
		return path_index == std::numeric_limits<SolnIndT>::max();
	}

	/// \brief Make the sentinel "no more work" task.
	static StartPointTask sentinel()
	{
		return StartPointTask{};  // default path_index == max
	}

	/// \cond PATH_RESULT_SERIALIZATION
	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & start_point;
	}
	/// \endcond
};


/**
\brief Result of executing one WHOLE path: start -> endgame boundary -> target.

The single result type for the speculative-full-path model.  A worker carries a path through both
the pre-endgame tracking and the endgame, so this bundles the boundary data (needed for the manager's
midpath/crossing check) together with the final solution and all endgame metadata.  The work item
that produces it is a StartPointTask (path index + rank 0's start point).
*/
template<typename ComplexT>
struct FullPathResult
{
	using SolnIndT = std::size_t;  ///< The path-index type.
	using RealT    = typename NumTraits<ComplexT>::Real;  ///< The real companion of the complex type.

	SolnIndT      path_index             = 0;  ///< The index of the path this result is for.

	// boundary (pre-endgame) data
	SuccessCode   pre_endgame_success_code    = SuccessCode::NeverStarted;  ///< Outcome of pre-endgame tracking to the boundary.
	Vec<ComplexT> boundary_point;  ///< The space point at the endgame boundary.
	RealT         boundary_stepsize      = RealT(0);  ///< The step size at the endgame boundary.
	unsigned      boundary_precision     = DoublePrecision();  ///< The precision at the endgame boundary.

	// endgame data
	SuccessCode   endgame_success_code        = SuccessCode::NeverStarted;  ///< Outcome of the endgame.
	Vec<ComplexT> final_solution;  ///< The final solution point.
	double        function_residual              = 0;  ///< Residual of the system at the final solution.
	double        condition_number               = 0;  ///< Condition-number estimate at the final solution.
	double        newton_residual                = 0;  ///< Final Newton residual.
	ComplexT      final_time_used;  ///< The time value the endgame finished at.
	double        accuracy_estimate              = 0;  ///< Estimated accuracy of the final solution.
	double        accuracy_estimate_user_coords  = 0;  ///< Estimated accuracy in the user's coordinates.
	unsigned      cycle_num                      = 0;  ///< The estimated cycle number at the endpoint.
	unsigned      precision_digits             = 0;   ///< Digits the endgame finished in (= solution point's precision).
	unsigned      accuracy_digits              = 0;   ///< Trustworthy digit count, from the convergence agreement.

	// precision metadata (spans the whole path)
	bool          precision_changed              = false;  ///< Whether precision changed during the path.
	ComplexT      time_of_first_prec_increase;  ///< The time of the first precision increase (if any).
	unsigned      max_precision_used             = 0;  ///< The highest precision used over the whole path.

	// wall-clock time to execute the whole path (pre-endgame + endgame), in seconds
	double        path_time_seconds              = 0;  ///< Wall-clock seconds to execute the whole path.

	/// \cond PATH_RESULT_SERIALIZATION
	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & pre_endgame_success_code;
		ar & boundary_point;
		ar & boundary_stepsize;
		ar & boundary_precision;
		ar & endgame_success_code;
		ar & final_solution;
		ar & function_residual;
		ar & condition_number;
		ar & newton_residual;
		ar & final_time_used;
		ar & accuracy_estimate;
		ar & accuracy_estimate_user_coords;
		ar & cycle_num;
		ar & precision_digits;
		ar & accuracy_digits;
		ar & precision_changed;
		ar & time_of_first_prec_increase;
		ar & max_precision_used;
		ar & path_time_seconds;
	}
	/// \endcond
};

namespace detail {

// Sentinel detection and factory for the work-item type.  The work item is a StartPointTask whose
// max path_index marks "no more work".  (The plain-size_t overloads remain for any other caller.)

/// \brief Query whether a plain index value is the sentinel.
inline bool is_sentinel(std::size_t v)
{
	return v == std::numeric_limits<std::size_t>::max();
}

/// \brief Make the sentinel index value.
inline std::size_t make_sentinel(std::size_t)
{
	return std::numeric_limits<std::size_t>::max();
}

/// \brief Query whether a work-item task is the sentinel.
template<typename ComplexT>
bool is_sentinel(StartPointTask<ComplexT> const& t)
{
	return t.is_sentinel();
}

/// \brief Make the sentinel work-item task.
template<typename ComplexT>
StartPointTask<ComplexT> make_sentinel(StartPointTask<ComplexT> const&)
{
	return StartPointTask<ComplexT>::sentinel();
}

} // namespace detail

} // namespace parallel
} // namespace bertini
