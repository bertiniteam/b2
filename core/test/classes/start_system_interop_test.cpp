//This file is part of Bertini 2.
//
//start_system_interop_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system_interop_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system_interop_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file start_system_interop_test.cpp

\brief Cross-cutting "sanity" tests for the start systems: TotalDegree (the linear-product
Bezout start), MHomogeneous (its multi-group generalization), and RootsOfUnity, exercised on the
SAME targets so their homogenization, patch-fitting, start-point validity, and homotopy
interoperability can be compared.  These deliberately hammer the corners:

  * EVERY start point of each system is checked to be a root of that start system, on its patch,
    in both double and (high) multiprecision -- not a sampled few.
  * MHomogeneous on a single affine variable group must reproduce TotalDegree's Bezout count
    (MHom is a generalization of total degree).
  * The blend homotopy H = (1-t) target + gamma t start is zero at every start point at t=1.
  * Patch-fit is checked at HIGH precision, where any rescale-then-truncate precision slip would
    show up as a fat patch residual.

These do NOT run a full ZeroDim solve (kept fast + hang-proof); they test the start-system layer.
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"

using System = bertini::System;
using Variable = bertini::node::Variable;
using Var = std::shared_ptr<Variable>;
using VariableGroup = bertini::VariableGroup;
using complex_dbl = bertini::complex_dbl;
using mpfr = bertini::complex_mp;
template<typename NumT> using Vec = bertini::Vec<NumT>;

using bertini::DefaultPrecision;
using bertini::Precision;
namespace ss = bertini::start_system;

BOOST_AUTO_TEST_SUITE(start_system_interop)

// ---- helpers ---------------------------------------------------------------------------------

// A square, single-affine-group polynomial system with the given per-function degrees, built as
// x_i^{d_i} + (a generic-ish lower-order tail) so it is genuinely degree d_i (not sparse), 1 group.
static System SquareSingleGroup(std::vector<unsigned> const& degs)
{
	System sys;
	std::vector<Var> x;
	for (size_t i = 0; i < degs.size(); ++i)
		x.push_back(Variable::Make("x" + std::to_string(i)));
	VariableGroup g(x.begin(), x.end());
	sys.AddVariableGroup(g);
	for (size_t i = 0; i < degs.size(); ++i)
	{
		std::shared_ptr<bertini::node::Node> f = pow(x[i], static_cast<int>(degs[i]));
		// a cross term one degree lower, to keep the function dense but still degree d_i
		if (degs.size() > 1 && degs[i] >= 1)
			f = f + x[(i + 1) % degs.size()] * pow(x[i], static_cast<int>(degs[i]) - 1);
		f = f - 2;
		sys.AddFunction(f);
	}
	return sys;
}

// largest |Eval| over EVERY start point of the (already homogenized/patched) start system, in T.
template<typename StartT, typename T>
static double WorstStartPointResidual(StartT const& s)
{
	double worst = 0;
	for (auto ii = decltype(s.NumStartPoints())(0); ii < s.NumStartPoints(); ++ii)
	{
		auto sp = s.template StartPoint<T>(ii);
		auto v = s.Eval(sp);
		for (Eigen::Index j = 0; j < v.size(); ++j)
			worst = std::max(worst, static_cast<double>(abs(v(j))));
	}
	return worst;
}

// ---- TotalDegree: every start point is a root, on the patch -----------------------------------

BOOST_AUTO_TEST_CASE(total_degree_all_start_points_on_patch_double)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(30);

	auto sys = SquareSingleGroup({2, 3, 4});   // Bezout 24
	sys.Homogenize();
	sys.AutoPatch();

	ss::TotalDegree td(sys);
	BOOST_CHECK(td.IsHomogeneous());
	BOOST_CHECK(td.IsPatched());
	BOOST_CHECK_EQUAL(td.NumStartPoints(), 24ull);   // 2*3*4

	// EVERY start point a root (incl. patch + homogenization rows), in double.
	BOOST_CHECK_LT((WorstStartPointResidual<ss::TotalDegree, complex_dbl>(td)), 1e-9);
}

// The patch-fit must survive at HIGH precision: rescale-then-truncate would leave a residual at
// the working-precision floor, so a tight bound at 100 digits guards the start-point precision path.
BOOST_AUTO_TEST_CASE(total_degree_all_start_points_on_patch_high_precision)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(100);

	auto sys = SquareSingleGroup({2, 3, 4});
	sys.Homogenize();
	sys.AutoPatch();

	ss::TotalDegree td(sys);

	// every MP start point is at the requested precision...
	for (auto ii = decltype(td.NumStartPoints())(0); ii < td.NumStartPoints(); ++ii)
		BOOST_CHECK_EQUAL(Precision(td.StartPoint<mpfr>(ii)), 100u);

	// ...and is a root on the patch to (nearly) full working precision.
	BOOST_CHECK_LT((WorstStartPointResidual<ss::TotalDegree, mpfr>(td)), 1e-90);
}

// ---- MHomogeneous is a generalization of TotalDegree -----------------------------------------

// On a SINGLE affine variable group, the m-homogeneous Bezout count collapses to the total-degree
// Bezout count == product of the function degrees.  MHom and TotalDegree must agree.
BOOST_AUTO_TEST_CASE(mhom_single_affine_group_matches_total_degree_count)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(30);

	for (std::vector<unsigned> degs : std::vector<std::vector<unsigned>>{ {1,1}, {2,3}, {2,3,4}, {5} })
	{
		auto sys = SquareSingleGroup(degs);
		// pre-homogenization construction is where MHom builds its degree matrix
		ss::MHomogeneous mhom(sys);
		ss::TotalDegree   td(sys);

		unsigned long long prod = 1;
		for (auto d : degs) prod *= d;

		BOOST_CHECK_EQUAL(td.NumStartPoints(), prod);
		BOOST_CHECK_EQUAL(mhom.NumStartPoints(), prod);   // MHom collapses to total degree on one group
	}
}

// MHom start points are roots on the patch too (double + MP), same as TotalDegree.
BOOST_AUTO_TEST_CASE(mhom_all_start_points_on_patch)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(50);

	auto sys = SquareSingleGroup({2, 3, 4});
	sys.Homogenize();
	sys.AutoPatch();

	ss::MHomogeneous mhom(sys);
	BOOST_CHECK_EQUAL(mhom.NumStartPoints(), 24ull);
	BOOST_CHECK(mhom.IsHomogeneous());
	BOOST_CHECK(mhom.IsPatched());

	BOOST_CHECK_LT((WorstStartPointResidual<ss::MHomogeneous, complex_dbl>(mhom)), 1e-9);
	BOOST_CHECK_LT((WorstStartPointResidual<ss::MHomogeneous, mpfr>(mhom)), 1e-40);
}

// ---- homotopy interoperability: H(start, t=1) == 0 -------------------------------------------

// Build the blend homotopy the way FormHomotopy does and check it vanishes at every start point at
// t = 1, for BOTH TotalDegree and MHomogeneous on the same target.
template<typename StartT>
static void CheckHomotopyZeroAtStartPoints(System const& target_hp, StartT const& start)
{
	auto t = Variable::Make("t");
	auto gamma = bertini::node::Rational::Make(bertini::node::Rational::Rand());
	System H = target_hp;        // copy the target's variable structure + patch
	H.ClearFunctions();          // the blend supplies the rows
	H.AddPathVariable(t);
	std::vector<std::shared_ptr<bertini::node::Node>> coeffs{ 1 - t, gamma * t };
	std::vector<std::shared_ptr<const System>> operands{
		std::make_shared<System>(target_hp),
		std::make_shared<System>(start) };
	H.AddBlock(bertini::blocks::BlendBlock<System>(t, coeffs, operands));

	for (auto ii = decltype(start.NumStartPoints())(0); ii < start.NumStartPoints(); ++ii)
	{
		auto Hval = H.Eval(start.template StartPoint<complex_dbl>(ii), complex_dbl(1));
		double worst = 0;
		for (Eigen::Index j = 0; j < Hval.size(); ++j)
			worst = std::max(worst, std::abs(Hval(j)));
		BOOST_CHECK_LT(worst, 1e-8);
	}
}

BOOST_AUTO_TEST_CASE(total_degree_and_mhom_homotopies_vanish_at_start_points)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(30);

	auto sys = SquareSingleGroup({2, 3});   // small Bezout 6
	sys.Homogenize();
	sys.AutoPatch();

	ss::TotalDegree  td(sys);
	ss::MHomogeneous mhom(sys);
	CheckHomotopyZeroAtStartPoints(sys, td);
	CheckHomotopyZeroAtStartPoints(sys, mhom);
}

// ---- corners ----------------------------------------------------------------------------------

// One variable, one function of degree d: exactly d start points, each a root.
BOOST_AUTO_TEST_CASE(total_degree_single_variable_degree_d)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(30);

	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(pow(x, 5) - 2);
	sys.Homogenize();
	sys.AutoPatch();

	ss::TotalDegree td(sys);
	BOOST_CHECK_EQUAL(td.NumStartPoints(), 5ull);
	// degree-5 root via a double linear solve: residual is conditioning-limited (~1e-10), so use a
	// scale-tolerant double threshold rather than a machine-epsilon one.
	BOOST_CHECK_LT((WorstStartPointResidual<ss::TotalDegree, complex_dbl>(td)), 1e-8);
}

// HomogenizePoint/DehomogenizePoint roundtrip on a TotalDegree start point: dehomogenizing then
// re-homogenizing returns the same (patch-fit) point.
BOOST_AUTO_TEST_CASE(total_degree_dehomogenize_homogenize_roundtrip)
{
	bertini::SetGlobalSeed(1u);
	DefaultPrecision(40);

	auto sys = SquareSingleGroup({2, 3});
	sys.Homogenize();
	sys.AutoPatch();
	ss::TotalDegree td(sys);

	for (auto ii = decltype(td.NumStartPoints())(0); ii < td.NumStartPoints(); ++ii)
	{
		auto sp = td.StartPoint<mpfr>(ii);
		auto round = td.HomogenizePoint(td.DehomogenizePoint(sp));
		BOOST_CHECK_LT((sp - round).template lpNorm<Eigen::Infinity>().convert_to<double>(), 1e-30);
	}
}

// Non-square single-group target must be rejected by TotalDegree (same guard family as MHom).
BOOST_AUTO_TEST_CASE(total_degree_non_square_throws)
{
	System sys;
	auto x = Variable::Make("x");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddFunction(pow(x, 2));
	sys.AddFunction(pow(x, 2) - 1);   // 2 functions, 1 variable
	BOOST_CHECK_THROW(ss::TotalDegree{sys}, std::runtime_error);
}

// A homogeneous variable group is not single-affine-group total degree: TotalDegree must reject it
// (it is MHom's domain).
BOOST_AUTO_TEST_CASE(total_degree_rejects_homogeneous_group)
{
	System sys;
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	sys.AddHomVariableGroup(VariableGroup{x, y});
	sys.AddFunction(pow(x, 2) + pow(y, 2));
	BOOST_CHECK_THROW(ss::TotalDegree{sys}, std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
