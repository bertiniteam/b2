//This file is part of Bertini 2.
//
//expand_to_function_tree_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//expand_to_function_tree_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with expand_to_function_tree_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file expand_to_function_tree_test.cpp

\brief Tests for System::ExpandToFunctionTree() -- the pure-function-tree twin of a
block-composed system.  The contract is that the expanded system evaluates and differentiates
identically to the block system, at every point.  This is the verification oracle used by the
AMP precision-escalation investigation to compare the block-composed homotopy's Jacobian against
an identical function-tree homotopy.
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/products_of_linears_block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"
#include "bertini2/function_tree.hpp"

BOOST_AUTO_TEST_SUITE(expand_to_function_tree)

using System        = bertini::System;
using Variable      = bertini::node::Variable;
using VariableGroup = bertini::VariableGroup;
using dbl           = bertini::dbl;
using mpfr_complex  = bertini::mpfr_complex;
template <typename T> using Vec = bertini::Vec<T>;
template <typename T> using Mat = bertini::Mat<T>;

namespace {

	// scale-relative comparison: |a-b| <= tol * (1 + |a|), entrywise.
	template <typename Derived>
	void CheckClose(Eigen::MatrixBase<Derived> const& a, Eigen::MatrixBase<Derived> const& b, double tol)
	{
		BOOST_REQUIRE_EQUAL(a.rows(), b.rows());
		BOOST_REQUIRE_EQUAL(a.cols(), b.cols());
		for (Eigen::Index i = 0; i < a.rows(); ++i)
			for (Eigen::Index j = 0; j < a.cols(); ++j)
			{
				double err   = static_cast<double>(abs(a(i, j) - b(i, j)));
				double scale = 1.0 + static_cast<double>(abs(a(i, j)));
				BOOST_CHECK_SMALL(err / scale, tol);
			}
	}

	Vec<dbl> RandomVecD(int n)
	{
		return Vec<dbl>::Random(n);
	}

	// compare block vs expanded system on eval + Jacobian at a point (autonomous overload)
	template <typename T>
	void CompareEvalJac(System const& blk, System const& exp, Vec<T> const& p, double tol)
	{
		CheckClose(blk.Eval(p),     exp.Eval(p),     tol);
		CheckClose(blk.Jacobian(p), exp.Jacobian(p), tol);
	}

	// ... and with a path-variable value
	template <typename T>
	void CompareEvalJac(System const& blk, System const& exp, Vec<T> const& p, T const& t, double tol)
	{
		CheckClose(blk.Eval(p, t),     exp.Eval(p, t),     tol);
		CheckClose(blk.Jacobian(p, t), exp.Jacobian(p, t), tol);
	}

} // anonymous namespace


// A plain polynomial system is already a function tree; expansion must be a numeric identity.
BOOST_AUTO_TEST_CASE(polynomial_block_roundtrip)
{
	using namespace bertini;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x * y - 1);
	sys.AddFunction(pow(x, 2) + y - 3);

	System expanded = sys.ExpandToFunctionTree();

	for (int trial = 0; trial < 8; ++trial)
		CompareEvalJac<dbl>(sys, expanded, RandomVecD(2), 1e-13);

	DefaultPrecision(50);
	sys.precision(50);
	expanded.precision(50);
	for (int trial = 0; trial < 4; ++trial)
		CompareEvalJac<mpfr_complex>(sys, expanded, RandomOfUnits<mpfr_complex>(2), 1e-40);
}


// A products-of-linears block: f_i = prod_r (linear form).  The expansion rebuilds those
// products as nodes; block and expanded must agree on eval and Jacobian.
BOOST_AUTO_TEST_CASE(products_of_linears_block_matches)
{
	using namespace bertini;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");

	System sys;
	sys.AddVariableGroup(VariableGroup{x, y});

	// f0 = (2x + 3y - 1)(x - y + 4);  f1 = (x + 5)( -y + 2)( x + y )   (3 factors)
	std::vector<Mat<mpfr_complex>> factors;
	{
		Mat<mpfr_complex> M0(2, 3); // rows = factors, cols = (x, y, const)
		M0(0,0) = 2; M0(0,1) = 3;  M0(0,2) = -1;
		M0(1,0) = 1; M0(1,1) = -1; M0(1,2) = 4;
		factors.push_back(M0);

		Mat<mpfr_complex> M1(3, 3);
		M1(0,0) = 1; M1(0,1) = 0;  M1(0,2) = 5;
		M1(1,0) = 0; M1(1,1) = -1; M1(1,2) = 2;
		M1(2,0) = 1; M1(2,1) = 1;  M1(2,2) = 0;
		factors.push_back(M1);
	}
	sys.AddBlock(blocks::ProductsOfLinearsBlock(2, factors));

	System expanded = sys.ExpandToFunctionTree();

	for (int trial = 0; trial < 8; ++trial)
		CompareEvalJac<dbl>(sys, expanded, RandomVecD(2), 1e-12);

	DefaultPrecision(50);
	sys.precision(50);
	expanded.precision(50);
	for (int trial = 0; trial < 4; ++trial)
		CompareEvalJac<mpfr_complex>(sys, expanded, RandomOfUnits<mpfr_complex>(2), 1e-40);
}


// A blend block H = (1-t)*A + gamma*t*B.  The expansion recurses into the operand systems and
// rebuilds the t-weighted sum; block and expanded must agree on eval and Jacobian (and the
// path-variable dependence).
BOOST_AUTO_TEST_CASE(blend_block_matches)
{
	using namespace bertini;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto t = Variable::Make("t");

	// operand A and operand B share the SAME variable nodes (x,y).
	System A;
	A.AddVariableGroup(VariableGroup{x, y});
	A.AddFunction(x * y - 1);
	A.AddFunction(x + y);

	System B;
	B.AddVariableGroup(VariableGroup{x, y});
	B.AddFunction(pow(x, 2) - 1);
	B.AddFunction(pow(y, 2) - 1);

	auto gamma = node::Rational::Make(node::Rational::Rand());

	System H;
	H.AddVariableGroup(VariableGroup{x, y});
	H.AddPathVariable(t);
	std::vector<std::shared_ptr<node::Node>> coeffs{ 1 - t, gamma * t };
	std::vector<std::shared_ptr<const System>> operands{
		std::make_shared<System>(A),
		std::make_shared<System>(B) };
	H.AddBlock(blocks::BlendBlock<System>(t, std::move(coeffs), std::move(operands)));

	System expanded = H.ExpandToFunctionTree();

	for (int trial = 0; trial < 8; ++trial)
	{
		dbl tau = dbl(0.3, -0.2) * dbl(trial + 1);
		CompareEvalJac<dbl>(H, expanded, RandomVecD(2), tau, 1e-12);
	}

	DefaultPrecision(50);
	H.precision(50);
	expanded.precision(50);
	for (int trial = 0; trial < 4; ++trial)
	{
		Vec<mpfr_complex> p = RandomOfUnits<mpfr_complex>(2);
		mpfr_complex tau = RandomOfUnits<mpfr_complex>(1)(0);
		CompareEvalJac<mpfr_complex>(H, expanded, p, tau, 1e-40);
	}
}


BOOST_AUTO_TEST_SUITE_END()
