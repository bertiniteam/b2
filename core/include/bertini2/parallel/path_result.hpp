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

constexpr int TAG_WORK_ITEM = 1;
constexpr int TAG_RESULT    = 2;
constexpr int TAG_CAPACITY  = 3;  // worker -> manager: int, max tasks in flight on that rank



/**
\brief Work item for the speculative-full-path model: a path index plus its start point.

Rank 0 computes start points authoritatively and ships them to workers, so a worker never derives
its own (that, with authoritative pi, is what keeps start points -- and hence whole tracks --
identical across serial and distributed runs).
*/
template<typename ComplexT>
struct StartPointTask
{
	using SolnIndT = std::size_t;

	SolnIndT      path_index = std::numeric_limits<SolnIndT>::max();  // max = sentinel
	Vec<ComplexT> start_point;

	bool is_sentinel() const
	{
		return path_index == std::numeric_limits<SolnIndT>::max();
	}

	static StartPointTask sentinel()
	{
		return StartPointTask{};  // default path_index == max
	}

	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & start_point;
	}
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
	using SolnIndT = std::size_t;
	using RealT    = typename NumTraits<ComplexT>::Real;

	SolnIndT      path_index             = 0;

	// boundary (pre-endgame) data
	SuccessCode   pre_endgame_success_code    = SuccessCode::NeverStarted;
	Vec<ComplexT> boundary_point;
	RealT         boundary_stepsize      = RealT(0);
	unsigned      boundary_precision     = DoublePrecision();

	// endgame data
	SuccessCode   endgame_success_code        = SuccessCode::NeverStarted;
	Vec<ComplexT> final_solution;
	double        function_residual              = 0;
	double        condition_number               = 0;
	double        newton_residual                = 0;
	ComplexT      final_time_used;
	double        accuracy_estimate              = 0;
	double        accuracy_estimate_user_coords  = 0;
	unsigned      cycle_num                      = 0;
	unsigned      precision_digits             = 0;   // digits the endgame finished in (= solution point's precision)
	unsigned      accuracy_digits              = 0;   // trustworthy digit count, from the convergence agreement

	// precision metadata (spans the whole path)
	bool          precision_changed              = false;
	ComplexT      time_of_first_prec_increase;
	unsigned      max_precision_used             = 0;

	// wall-clock time to execute the whole path (pre-endgame + endgame), in seconds
	double        path_time_seconds              = 0;

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
};

namespace detail {

// Sentinel detection and factory for the work-item type.  The work item is a StartPointTask whose
// max path_index marks "no more work".  (The plain-size_t overloads remain for any other caller.)

inline bool is_sentinel(std::size_t v)
{
	return v == std::numeric_limits<std::size_t>::max();
}

inline std::size_t make_sentinel(std::size_t)
{
	return std::numeric_limits<std::size_t>::max();
}

template<typename ComplexT>
bool is_sentinel(StartPointTask<ComplexT> const& t)
{
	return t.is_sentinel();
}

template<typename ComplexT>
StartPointTask<ComplexT> make_sentinel(StartPointTask<ComplexT> const&)
{
	return StartPointTask<ComplexT>::sentinel();
}

} // namespace detail

} // namespace parallel
} // namespace bertini
