//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/newton_refine.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/newton_refine.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/newton_refine.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// silviana amethyst, university of wisconsin eau claire

/**
\file bertini2/nag_algorithms/newton_refine.hpp

\brief Standalone Newton refinement of a point against a System -- no tracker required.

The sharpening primitive: given a square system and an approximate solution, iterate
Newton's method until consecutive approximations agree to a requested tolerance.  The
system may be autonomous (no path variable -- the ordinary case for deflated critical
point systems) or carry a path variable with a caller-supplied time.

The whole point of the standalone (versus the tracker-owned \c Refine) is refinement on
DEFLATED systems at singular points: deflation restores quadratic convergence exactly
where the tracker's own system cannot converge, and a deflated system is a plain
autonomous System that no tracker is configured around.  Overdetermined systems are
refused with an instructive message -- square them by randomization first, as the
tracking layer does.

Evaluation runs at the precision of the supplied system and point; lift both before
calling to sharpen beyond their current precision.
*/

#pragma once

#include "bertini2/system/system.hpp"
#include "bertini2/trackers/config.hpp"
#include "bertini2/linalg/lu_solver.hpp"

namespace bertini {
namespace algorithm {

/**
\brief The outcome of a standalone Newton refinement.

\tparam ComplexT the complex number type the refinement ran in.
*/
template <typename ComplexT>
struct NewtonRefineResult
{
	tracking::SuccessCode code;  ///< Success when the step norm reached the tolerance; FailedToConverge or MatrixSolveFailure otherwise.
	Vec<ComplexT> point;         ///< The refined point (the best iterate reached, even on failure).
	NumErrorT achieved;          ///< Infinity norm of the last Newton step -- the consecutive-approximation agreement actually achieved.
	unsigned iterations;         ///< Number of Newton iterations taken.
};

/**
\brief Newton-refine a point against a square System, without a tracker.

Iterates full Newton steps until the infinity norm of the step falls at or below
\p tolerance, or \p max_iterations steps have been taken.  For a system with a path
variable, evaluation is at time \p time; for an autonomous system \p time is ignored.

\tparam ComplexT the complex number type to iterate in.

\param S the system to refine against.  Must be square: \c NumTotalFunctions()
	(user functions plus patches) equal to \c NumVariables().  Deflated systems are
	the intended customers -- square them by randomization if overdetermined.
\param start the approximate solution to refine.
\param tolerance stop when the infinity norm of a Newton step is at or below this.
\param max_iterations refuse to iterate more than this many times.
\param time the path-variable value for non-autonomous systems; ignored otherwise.

\return a NewtonRefineResult carrying the refined point, the achieved step norm,
	the iteration count, and the SuccessCode.

\throws std::runtime_error if the system is not square, with a message saying how
	to square it.
*/
template <typename ComplexT>
NewtonRefineResult<ComplexT> NewtonRefine(System const& S,
                                          Vec<ComplexT> const& start,
                                          NumErrorT tolerance,
                                          unsigned max_iterations,
                                          ComplexT time = ComplexT(0))
{
	const auto n_funcs = S.NumTotalFunctions();
	const auto n_vars = S.NumVariables();
	if (n_funcs != n_vars)
	{
		std::stringstream ss;
		ss << "NewtonRefine requires a SQUARE system, but this one has "
		   << n_funcs << " total functions (user functions plus patches) over "
		   << n_vars << " variables.  Square an overdetermined system by "
		   << "randomization (multiply by a generic full-rank matrix) before refining.";
		throw std::runtime_error(ss.str());
	}
	if (start.size() != static_cast<Eigen::Index>(n_vars))
	{
		std::stringstream ss;
		ss << "NewtonRefine: start point has " << start.size()
		   << " coordinates but the system has " << n_vars << " variables.";
		throw std::runtime_error(ss.str());
	}

	NewtonRefineResult<ComplexT> result{tracking::SuccessCode::FailedToConverge,
	                                    start, static_cast<NumErrorT>(-1), 0};

	Vec<ComplexT> f(n_funcs);
	Mat<ComplexT> J(n_funcs, n_vars);
	Vec<ComplexT> step(n_vars);
	linalg::PartialPivLU<ComplexT> lu;

	for (unsigned it = 0; it < max_iterations; ++it)
	{
		if (S.HavePathVariable())
			S.SetAndReset(result.point, time);
		else
			S.SetAndReset(result.point);
		S.EvalInPlace(f);
		S.JacobianInPlace(J);

		if (lu.Factor(J) != MatrixSuccessCode::Success)
		{
			result.code = tracking::SuccessCode::MatrixSolveFailure;
			return result;
		}
		lu.Solve(f, step);          // step = +J^{-1} f = -(Newton step)
		result.point -= step;
		result.iterations = it + 1;
		result.achieved = static_cast<NumErrorT>(
			step.template lpNorm<Eigen::Infinity>());

		if (result.achieved <= tolerance)
		{
			result.code = tracking::SuccessCode::Success;
			return result;
		}
	}
	return result;
}

} // namespace algorithm
} // namespace bertini
