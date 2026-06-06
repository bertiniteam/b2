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

\brief Task and result structs used in the manager-worker MPI protocol for ZeroDim.

Phase 1 (before endgame):
  task   = SolnIndT path index
  result = PathBeforeEGResult<ComplexT>

Phase 2 (during endgame):
  task   = Phase2Task<ComplexT>  (includes boundary point so workers are self-contained)
  result = PathDuringEGResult<ComplexT>

All structs carry Boost.Serialization support so Boost.MPI can transmit them
via the existing serialization for Eigen vectors and arbitrary-precision types.
*/

#pragma once

#include <deque>
#include <limits>

#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"   // defines Vec<T> and sets up Eigen serialization plugin
#include "bertini2/mpfr_extensions.hpp"    // Boost.Serialization for mpfr_float, mpc_complex
#include "bertini2/common/config.hpp"      // SuccessCode (needs Vec<T> and <deque> first)

#include <boost/serialization/vector.hpp>

namespace bertini {
namespace parallel {

constexpr int TAG_WORK_ITEM = 1;
constexpr int TAG_RESULT    = 2;


/**
\brief Result of tracking a single path from start time to the endgame boundary.

Sent from worker to manager after TrackSinglePathBeforeEG(idx).
*/
template<typename ComplexT>
struct PathBeforeEGResult
{
	using SolnIndT = std::size_t;
	using RealT    = typename NumTraits<ComplexT>::Real;

	SolnIndT      path_index             = 0;
	SuccessCode   pre_endgame_success    = SuccessCode::NeverStarted;
	Vec<ComplexT> boundary_point;
	RealT         boundary_stepsize      = RealT(0);
	bool          precision_changed      = false;
	ComplexT      time_of_first_prec_increase;
	unsigned      max_precision_used     = 0;

	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & pre_endgame_success;
		ar & boundary_point;
		ar & boundary_stepsize;
		ar & precision_changed;
		ar & time_of_first_prec_increase;
		ar & max_precision_used;
	}
};


/**
\brief Task for Phase 2 (during-endgame) tracking.

Carries the path index and the boundary point/stepsize so each worker is
self-contained: workers only tracked a subset of paths in Phase 1 and may
not have the boundary data for paths assigned to them in Phase 2.
*/
template<typename ComplexT>
struct Phase2Task
{
	using SolnIndT = std::size_t;
	using RealT    = typename NumTraits<ComplexT>::Real;

	SolnIndT      path_index        = std::numeric_limits<SolnIndT>::max();  // max = sentinel
	Vec<ComplexT> boundary_point;
	RealT         boundary_stepsize = RealT(0);

	bool is_sentinel() const
	{
		return path_index == std::numeric_limits<SolnIndT>::max();
	}

	static Phase2Task sentinel()
	{
		return Phase2Task{};  // default path_index == max
	}

	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & boundary_point;
		ar & boundary_stepsize;
	}
};


/**
\brief Result of tracking a single path through the endgame.

Sent from worker to manager after TrackSinglePathDuringEG.
Carries the final solution and all endgame metadata fields that would
normally be written into SolutionMetaData by the tracking code.
*/
template<typename ComplexT>
struct PathDuringEGResult
{
	using SolnIndT = std::size_t;

	SolnIndT      path_index                     = 0;
	Vec<ComplexT> final_solution;

	SuccessCode   endgame_success                = SuccessCode::NeverStarted;
	double        function_residual              = 0;
	double        condition_number               = 0;
	double        newton_residual                = 0;
	ComplexT      final_time_used;
	double        accuracy_estimate              = 0;
	double        accuracy_estimate_user_coords  = 0;
	unsigned      cycle_num                      = 0;
	bool          precision_changed              = false;
	ComplexT      time_of_first_prec_increase;
	unsigned      max_precision_used             = 0;

	template<class Archive>
	void serialize(Archive& ar, unsigned const)
	{
		ar & path_index;
		ar & final_solution;
		ar & endgame_success;
		ar & function_residual;
		ar & condition_number;
		ar & newton_residual;
		ar & final_time_used;
		ar & accuracy_estimate;
		ar & accuracy_estimate_user_coords;
		ar & cycle_num;
		ar & precision_changed;
		ar & time_of_first_prec_increase;
		ar & max_precision_used;
	}
};

namespace detail {

// Sentinel detection and factory — overloaded for each task type so manager
// and worker can use the same logic without knowing the concrete type.

inline bool is_sentinel(std::size_t v)
{
	return v == std::numeric_limits<std::size_t>::max();
}

inline std::size_t make_sentinel(std::size_t)
{
	return std::numeric_limits<std::size_t>::max();
}

template<typename ComplexT>
bool is_sentinel(Phase2Task<ComplexT> const& t)
{
	return t.is_sentinel();
}

template<typename ComplexT>
Phase2Task<ComplexT> make_sentinel(Phase2Task<ComplexT> const&)
{
	return Phase2Task<ComplexT>::sentinel();
}

} // namespace detail

} // namespace parallel
} // namespace bertini
