//This file is part of Bertini 2.
//
//python/linalg_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/linalg_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/linalg_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin - Eau Claire
//  Summer 2026
//
//
//  python/linalg_export.cpp:  Dense linear algebra for the multiprecision types.
//
//  The trick here is that eigenpy's matrix decompositions (LU, QR, ...) are exposed by
//  header-only def_visitor templates parameterized on the matrix type -- e.g.
//  eigenpy::PartialPivLUSolverVisitor<MatrixType>.  A stock `import eigenpy` only baked in
//  instantiations for the standard scalars (double / complex128) and cannot make an mp one at
//  runtime.  But THIS translation unit is compiled with knowledge of both eigenpy's headers and
//  bertini's mp scalars, so it can instantiate those same visitors on Mat<complex_mp> /
//  Mat<real_mp>.  That reuses eigenpy's own binding code + Eigen's algorithm -- we do not
//  reimplement LU -- and the resulting classes live under bertini.linalg (they only register when
//  _pybertini is imported, so bertini is their honest home).

#include "python_common.hpp"

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"

#include <eigenpy/decompositions/PartialPivLU.hpp>

namespace bertini{
	namespace python{

		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

		namespace {

			template<typename T>
			Vec<T> SolveLinearSystem(Mat<T> const& A, Vec<T> const& b)
			{
				return Vec<T>(A.partialPivLu().solve(b));
			}

		} // anonymous namespace


		void ExportLinalg()
		{
			using namespace boost::python;

			// Build the bertini.linalg submodule (same idiom as bertini.multiprec / bertini.function_tree).
			scope current_scope;
			std::string submodule_name(extract<const char*>(current_scope.attr("__name__")));
			submodule_name.append(".linalg");
			object linalg_module(borrowed(PyImport_AddModule(submodule_name.c_str())));
			current_scope.attr("linalg") = linalg_module;

			scope linalg_scope = linalg_module;
			linalg_scope.attr("__doc__") =
				"Dense linear algebra for bertini's multiprecision types (real_mp / complex_mp), "
				"at full multiprecision.  Backed by eigenpy's own Eigen decomposition wrappers "
				"instantiated on the mp scalars -- the LU that a stock `import eigenpy` cannot do "
				"on these custom types.";

			// eigenpy's own PartialPivLU visitor, instantiated on the mp matrix types.
			eigenpy::PartialPivLUSolverVisitor<Mat<complex_mp>>::expose("PartialPivLU");
			eigenpy::PartialPivLUSolverVisitor<Mat<real_mp>>::expose("PartialPivLUReal");

			// One-shot convenience: solve a square system A x = b at mp precision, partial-pivot LU.
			def("solve",
				+[](Mat<complex_mp> const& A, Vec<complex_mp> const& b){ return SolveLinearSystem<complex_mp>(A, b); },
				(arg("A"), arg("b")),
				"Solve the square linear system A x = b at multiprecision (complex_mp), via partial-pivot LU.  Returns x.");
			def("solve",
				+[](Mat<real_mp> const& A, Vec<real_mp> const& b){ return SolveLinearSystem<real_mp>(A, b); },
				(arg("A"), arg("b")),
				"Solve the square linear system A x = b at multiprecision (real_mp), via partial-pivot LU.  Returns x.");
		}

	} // namespace python
} // namespace bertini
