//This file is part of Bertini 2.
//
//describe.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//describe.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with describe.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/system/blocks/describe.hpp

\brief Small shared helpers for the per-block `Describe` methods -- the human-facing printing of a
block-composed System.  Each block prints the rows it owns, labelled by a running global index
`f_k`.  Two levels: **terse** (default) uses placeholder symbols for matrices/coefficients so only
the structure shows; **verbose** reveals the actual numbers.  This output is for reading, not for
re-parsing.
*/

#pragma once

#include <ostream>
#include <string>

#include "bertini2/num_traits.hpp"
#include "bertini2/function_tree.hpp"

namespace bertini {
namespace blocks {
namespace describe_detail {

/// How many function rows a block lists in *terse* (non-verbose) describe before truncating with a
/// `... (k more ...)` line.  Keeps a large system (e.g. a slice with many forms) from flooding the
/// terminal now that the terse output carries actual coefficients; verbose prints every row.
inline constexpr size_t kTerseRowCap = 10;

/// Format one multiprecision real part.  `sig > 0` renders that many significant digits (the terse
/// default -- short and legible); `sig <= 0` renders full precision (verbose / round-trippable).
template <typename R>
inline std::string FmtReal(R const& r, int sig)
{
	return (sig > 0) ? r.str(sig) : r.str();
}

/// Print one (complex_mp) coefficient compactly: a real prints as its real part, a pure imaginary
/// as `b*i`, otherwise as `(a+b*i)`.  `sig` controls the digit count (see FmtReal); the default of
/// 0 keeps full precision for any caller that does not opt into the short form.
inline void PrintCoeff(std::ostream& out, complex_mp const& c, int sig = 0)
{
	const bool re0 = (c.real() == 0), im0 = (c.imag() == 0);
	if (im0)            out << FmtReal(c.real(), sig);
	else if (re0)       out << FmtReal(c.imag(), sig) << "*i";
	else                out << "(" << FmtReal(c.real(), sig) << "+" << FmtReal(c.imag(), sig) << "*i)";
}

/// The significant-digit count for coefficient printing: 4 for terse (short and legible), the
/// current working precision for verbose (every digit the system actually carries -- the master
/// coefficient matrix is stored at a much higher precision than that, so printing its full string
/// would dump thousands of digits).
inline int CoeffSig(bool verbose) { return verbose ? static_cast<int>(DefaultPrecision()) : 4; }

/// `f_k` for a single row, or `f_a..f_b` for a contiguous range of `n` rows starting at `row`.
inline void PrintRowLabel(std::ostream& out, size_t row, size_t n)
{
	out << "f_" << row;
	if (n > 1)
		out << "..f_" << (row + n - 1);
}

/// The augmenting variable list `[x, y, 1]` (affine) or `[h, x, y]` (homogeneous: no trailing 1).
inline void PrintAugmentedVars(std::ostream& out, VariableGroup const& vars, size_t num_vars, bool homogeneous)
{
	out << "[";
	for (size_t c = 0; c < num_vars && c < vars.size(); ++c)
		out << (c ? ", " : "") << *vars[c];
	if (!homogeneous)
		out << ", 1";
	out << "]";
}

/// One affine linear form, row r of an augmented coefficient matrix M (num_vars+1 cols affine, or
/// num_vars cols homogeneous).  Verbose prints the actual sum `2*x + 1*y - 1`; otherwise nothing
/// (the caller prints the placeholder `c.[...]`).
inline void PrintLinearFormVerbose(std::ostream& out, Mat<complex_mp> const& M, Eigen::Index r,
                                   VariableGroup const& vars, size_t num_vars, bool homogeneous)
{
	bool first = true;
	const size_t ncol = homogeneous ? num_vars : num_vars + 1;
	for (size_t c = 0; c < ncol; ++c)
	{
		complex_mp const& coeff = M(r, static_cast<Eigen::Index>(c));
		if (coeff.real() == 0 && coeff.imag() == 0)
			continue;
		if (!first) out << " + ";
		first = false;
		out << "(";
		PrintCoeff(out, coeff);
		out << ")";
		if (c < num_vars && c < vars.size())     // a variable column (the last affine column is the constant)
			out << "*" << *vars[c];
	}
	if (first)
		out << "0";
}

} // namespace describe_detail
} // namespace blocks
} // namespace bertini
