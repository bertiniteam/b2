//This file is part of Bertini 2.
//
//python/newton_refine_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/newton_refine_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/newton_refine_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// silviana amethyst, university of wisconsin eau claire

#include "newton_refine_export.hpp"

#include "bertini2/nag_algorithms/newton_refine.hpp"
#include "bertini2/system/system.hpp"

namespace bertini {
namespace python {

namespace {

// ADR-0001: the scalar time is taken BY VALUE (never const& adjacent to an
// eigenpy-converted vector), and the refined point is RETURNED rather than
// written through a writable Ref.
boost::python::tuple NewtonRefineAtTime(bertini::System const& sys,
                                        Vec<complex_mp> const& start,
                                        double tolerance,
                                        unsigned max_iterations,
                                        complex_mp time)
{
	auto result = bertini::algorithm::NewtonRefine(sys, start, tolerance,
	                                               max_iterations, time);
	return boost::python::make_tuple(result.point, result.code,
	                                 result.achieved, result.iterations);
}

boost::python::tuple NewtonRefineAutonomous(bertini::System const& sys,
                                            Vec<complex_mp> const& start,
                                            double tolerance,
                                            unsigned max_iterations)
{
	return NewtonRefineAtTime(sys, start, tolerance, max_iterations,
	                          complex_mp(0));
}

} // namespace


void ExportNewtonRefine()
{
	boost::python::def(
		"newton_refine", &NewtonRefineAutonomous,
		"newton_refine(system, start, tolerance, max_iterations) -> (point, code, achieved, iterations)\n\n"
		"Standalone Newton refinement of a point against a SQUARE system -- no tracker.\n"
		"Iterates full Newton steps until consecutive approximations agree to `tolerance`\n"
		"in the infinity norm, at the precision of the supplied system and point.\n\n"
		"The intended customers are DEFLATED systems at singular points: deflation\n"
		"restores quadratic convergence exactly where the plain system cannot converge.\n"
		"Overdetermined systems raise with instructions to square by randomization.\n\n"
		"Returns (refined point, SuccessCode, achieved step norm, iterations taken).");

	boost::python::def(
		"newton_refine", &NewtonRefineAtTime,
		"newton_refine(system, start, tolerance, max_iterations, time) -> (point, code, achieved, iterations)\n\n"
		"As newton_refine/4, for a system with a path variable, evaluated at `time`.");
}

} // namespace python
} // namespace bertini
