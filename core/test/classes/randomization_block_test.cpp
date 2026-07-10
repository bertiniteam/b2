//This file is part of Bertini 2.
//
//randomization_block_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//randomization_block_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file randomization_block_test.cpp

\brief Unit tests for RandomizationBlock and System::Randomize.

The workhorse oracle is ExpandToFunctionTree(): a randomized system evaluated through its block
must agree, value-for-value and derivative-for-derivative, with the same system expanded to plain
function-tree nodes -- in both complex_dbl and complex_mp, before AND after homogenization (so the
homogenizing-variable power deficits are exercised), for single- and multi-projective systems.
That cross-check is deterministic on every platform, unlike a heap-dirtiness-dependent bug.
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/randomization_block.hpp"

BOOST_AUTO_TEST_SUITE(randomization_block_suite)

using namespace bertini;
using bertini::node::Variable;
using bertini::blocks::RandomizationBlock;
using Randomization = RandomizationBlock<System>;

// The contract is checkable here (System is complete), even though block.hpp cannot static_assert
// it (System is only forward-declared there).
static_assert(bertini::blocks::is_block_v<Randomization>,
              "RandomizationBlock<System> must satisfy the block contract");


// ---- helpers -------------------------------------------------------------------------------

// An overdetermined single-affine-group system of UNEQUAL degrees (so randomization needs the
// h-power deficits once homogenized): f0 = x^2 + y^2 - 1 (deg 2), f1 = x - y (deg 1),
// f2 = 2x^2 - 1 (deg 2).  Common zeros are the two points (+/- 1/sqrt2, +/- 1/sqrt2) with x = y.
static System OverdeterminedSingleGroup()
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System s;
	s.AddVariableGroup(VariableGroup{x, y});
	s.AddFunction(x * x + y * y - node::Integer::Make(1));
	s.AddFunction(x - y);
	s.AddFunction(node::Integer::Make(2) * x * x - node::Integer::Make(1));
	return s;
}

// An overdetermined two-affine-group system of UNEQUAL multidegrees: with groups {x}, {y},
// f0 = x*y (md (1,1)), f1 = x*x*y (md (2,1)), f2 = x (md (1,0)), f3 = y (md (0,1)).
static System OverdeterminedMultiGroup()
{
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	System s;
	s.AddVariableGroup(VariableGroup{x});
	s.AddVariableGroup(VariableGroup{y});
	s.AddFunction(x * y);
	s.AddFunction(x * x * y);
	s.AddFunction(x);
	s.AddFunction(y);
	return s;
}

// Compare a block-composed system against its function-tree expansion at a handful of points, in
// both numeric types.  Covers values AND the Jacobian.
static void AgreesWithExpansion(System const& sys)
{
	System twin = sys.ExpandToFunctionTree();

	BOOST_REQUIRE_EQUAL(sys.NumNaturalFunctions(), twin.NumNaturalFunctions());
	BOOST_REQUIRE_EQUAL(sys.NumVariables(), twin.NumVariables());
	// the randomized rows' degrees must survive the expansion
	BOOST_CHECK(sys.Degrees() == twin.Degrees());

	const Eigen::Index nv = static_cast<Eigen::Index>(sys.NumVariables());

	// a few non-degenerate complex points
	std::vector<Vec<complex_dbl>> pts;
	{
		Vec<complex_dbl> p(nv);
		for (Eigen::Index k = 0; k < nv; ++k) p(k) = complex_dbl(0.3 + 0.17 * static_cast<double>(k), 0.5 - 0.11 * static_cast<double>(k));
		pts.push_back(p);
		Vec<complex_dbl> q(nv);
		for (Eigen::Index k = 0; k < nv; ++k) q(k) = complex_dbl(-0.7 + 0.05 * static_cast<double>(k), 0.9 + 0.03 * static_cast<double>(k));
		pts.push_back(q);
	}

	for (auto const& p : pts)
	{
		// --- complex_dbl ---
		Vec<complex_dbl> ea = sys.Eval(p);
		Vec<complex_dbl> eb = twin.Eval(p);
		BOOST_REQUIRE_EQUAL(ea.size(), eb.size());
		for (Eigen::Index i = 0; i < ea.size(); ++i)
			BOOST_CHECK(std::abs(ea(i) - eb(i)) < 1e-10);

		Mat<complex_dbl> ja = sys.Jacobian(p);
		Mat<complex_dbl> jb = twin.Jacobian(p);
		BOOST_REQUIRE_EQUAL(ja.rows(), jb.rows());
		BOOST_REQUIRE_EQUAL(ja.cols(), jb.cols());
		for (Eigen::Index i = 0; i < ja.rows(); ++i)
			for (Eigen::Index j = 0; j < ja.cols(); ++j)
				BOOST_CHECK(std::abs(ja(i, j) - jb(i, j)) < 1e-10);

		// --- complex_mp ---
		DefaultPrecision(40);
		Vec<complex_mp> pm(nv);
		for (Eigen::Index k = 0; k < nv; ++k)
			pm(k) = complex_mp(p(k).real(), p(k).imag());

		Vec<complex_mp> ema = sys.Eval(pm);
		Vec<complex_mp> emb = twin.Eval(pm);
		for (Eigen::Index i = 0; i < ema.size(); ++i)
			BOOST_CHECK(abs(ema(i) - emb(i)) < real_mp("1e-30"));

		Mat<complex_mp> jma = sys.Jacobian(pm);
		Mat<complex_mp> jmb = twin.Jacobian(pm);
		for (Eigen::Index i = 0; i < jma.rows(); ++i)
			for (Eigen::Index j = 0; j < jma.cols(); ++j)
				BOOST_CHECK(abs(jma(i, j) - jmb(i, j)) < real_mp("1e-30"));
	}
}


// ---- shape / construction ------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(randomize_squares_the_system)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();
	BOOST_REQUIRE_EQUAL(s.NumNaturalFunctions(), 3u);

	System r = s.Randomize();

	BOOST_CHECK_EQUAL(r.NumNaturalFunctions(), 2u);   // n = NumVariables - NumHomVariableGroups = 2
	BOOST_CHECK_EQUAL(r.NumVariables(), 2u);
	// the original system is NOT mutated
	BOOST_CHECK_EQUAL(s.NumNaturalFunctions(), 3u);
}

BOOST_AUTO_TEST_CASE(randomize_does_not_mutate_original_ordering)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();
	auto degs_before = s.Degrees();      // (2,1,2) in author order
	s.Randomize();
	auto degs_after = s.Degrees();
	BOOST_CHECK(degs_before == degs_after);  // the descending sort happened on the copy, not on s
}

BOOST_AUTO_TEST_CASE(optimal_path_count_is_product_of_top_n_degrees)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();   // degrees 2,1,2 -> top two are 2,2
	System r = s.Randomize();
	auto d = r.Degrees();
	BOOST_REQUIRE_EQUAL(d.size(), 2u);
	int product = 1;
	for (int e : d) product *= e;
	BOOST_CHECK_EQUAL(product, 4);            // 2*2, not 2^3 nor (maxdeg)^n with the low one inflated
}

BOOST_AUTO_TEST_CASE(matrix_getter_round_trips_and_is_I_C)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();
	System r = s.Randomize();

	auto const& block = std::get<Randomization>(r.Blocks().front());
	auto const& R = block.RandomizationMatrix();
	BOOST_REQUIRE_EQUAL(R.rows(), 2);
	BOOST_REQUIRE_EQUAL(R.cols(), 3);
	// leading 2x2 is the identity (the [I | C] structure after the descending sort)
	BOOST_CHECK(R(0, 0) == complex_mp(1) && R(1, 1) == complex_mp(1));
	BOOST_CHECK(R(0, 1) == complex_mp(0) && R(1, 0) == complex_mp(0));
}

BOOST_AUTO_TEST_CASE(underdetermined_throws)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	System s;
	s.AddVariableGroup(VariableGroup{x, y, z});
	s.AddFunction(x + y + z);            // 1 function, 3 variables
	BOOST_CHECK_THROW(s.Randomize(), std::runtime_error);
}


// ---- randomize to a chosen codimension -----------------------------------------------------

BOOST_AUTO_TEST_CASE(randomize_to_codimension_one)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();   // 3 functions, degrees (2,1,2)
	BOOST_REQUIRE_EQUAL(s.NumNaturalFunctions(), 3u);

	System r = s.Randomize(1);

	BOOST_CHECK_EQUAL(r.NumNaturalFunctions(), 1u);   // reduced to a single function
	auto d = r.Degrees();
	BOOST_REQUIRE_EQUAL(d.size(), 1u);
	BOOST_CHECK_EQUAL(d[0], 2);                        // leading (top-degree) function after the sort
	// the original system is NOT mutated
	BOOST_CHECK_EQUAL(s.NumNaturalFunctions(), 3u);
}

BOOST_AUTO_TEST_CASE(randomize_to_codimension_matches_square)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();           // n = NumVariables - NumHomVariableGroups = 2
	System r = s.Randomize(2);                          // explicit codimension == the affine dimension
	BOOST_CHECK_EQUAL(r.NumNaturalFunctions(), 2u);    // same shape as the no-arg square Randomize()
	BOOST_CHECK_EQUAL(s.Randomize().NumNaturalFunctions(), 2u);
}

BOOST_AUTO_TEST_CASE(randomize_codimension_zero_throws)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();
	BOOST_CHECK_THROW(s.Randomize(0), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(randomize_codimension_negative_throws)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();
	BOOST_CHECK_THROW(s.Randomize(-1), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(randomize_codimension_ge_num_functions_throws)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();            // 3 natural functions
	BOOST_CHECK_THROW(s.Randomize(3), std::runtime_error);   // == N: not a reduction
	BOOST_CHECK_THROW(s.Randomize(4), std::runtime_error);   // >  N
}

BOOST_AUTO_TEST_CASE(codimension_row_matches_expansion)
{
	DefaultPrecision(40);
	System r = OverdeterminedSingleGroup().Randomize(1);   // a single randomized function
	AgreesWithExpansion(r);
	r.Homogenize();                                        // exercises the h-power deficit
	AgreesWithExpansion(r);
}


// ---- evaluation correctness via the expansion oracle ---------------------------------------

BOOST_AUTO_TEST_CASE(affine_single_group_matches_expansion)
{
	DefaultPrecision(40);
	System r = OverdeterminedSingleGroup().Randomize();   // affine, not yet homogenized
	AgreesWithExpansion(r);
}

BOOST_AUTO_TEST_CASE(homogenized_single_group_matches_expansion)
{
	DefaultPrecision(40);
	System r = OverdeterminedSingleGroup().Randomize();
	r.Homogenize();                                       // exercises the h-power deficit (deg-1 fn padded by h)
	AgreesWithExpansion(r);
}

BOOST_AUTO_TEST_CASE(affine_multi_group_matches_expansion)
{
	DefaultPrecision(40);
	System r = OverdeterminedMultiGroup().Randomize();
	AgreesWithExpansion(r);
}

BOOST_AUTO_TEST_CASE(homogenized_multi_group_matches_expansion)
{
	DefaultPrecision(40);
	System r = OverdeterminedMultiGroup().Randomize();
	r.Homogenize();                                       // per-group hom-var power deficits
	AgreesWithExpansion(r);
}


// ---- user-supplied matrix ------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(user_supplied_matrix_reproduces_combination)
{
	DefaultPrecision(40);
	System s = OverdeterminedSingleGroup();    // f0=x^2+y^2-1, f1=x-y, f2=2x^2-1 (order preserved)

	// g0 = 1*f0 + 0*f1 + 0*f2 ; g1 = 3*f0 + 0*f1 + 5*f2  (kept in author order: f0 has max degree)
	Mat<complex_mp> R(2, 3);
	R << complex_mp(1), complex_mp(0), complex_mp(0),
	     complex_mp(3), complex_mp(0), complex_mp(5);
	System r = s.Randomize(R);

	BOOST_CHECK_EQUAL(r.NumNaturalFunctions(), 2u);

	// at (x,y) = (2,3):  f0 = 4+9-1 = 12, f2 = 8-1 = 7.  g0 = 12, g1 = 3*12 + 5*7 = 71.
	Vec<complex_dbl> p(2); p << complex_dbl(2), complex_dbl(3);
	Vec<complex_dbl> g = r.Eval(p);
	BOOST_CHECK(std::abs(g(0) - complex_dbl(12)) < 1e-10);
	BOOST_CHECK(std::abs(g(1) - complex_dbl(71)) < 1e-10);

	// and it agrees with its own expansion, affine and homogenized
	AgreesWithExpansion(r);
	r.Homogenize();
	AgreesWithExpansion(r);
}

BOOST_AUTO_TEST_CASE(user_supplied_matrix_wrong_columns_throws)
{
	DefaultPrecision(30);
	System s = OverdeterminedSingleGroup();   // 3 natural functions
	Mat<complex_mp> R(2, 2);                // 2 columns != 3
	R << complex_mp(1), complex_mp(0), complex_mp(0), complex_mp(1);
	BOOST_CHECK_THROW(s.Randomize(R), std::runtime_error);
}


// ---- ADR-0021: a block fully defines its own rows into a dirty buffer ----------------------

BOOST_AUTO_TEST_CASE(eval_and_jacobian_fully_define_rows_in_a_dirty_buffer)
{
	DefaultPrecision(40);
	System r = OverdeterminedSingleGroup().Randomize();
	r.Homogenize();

	auto const& block = std::get<Randomization>(r.Blocks().front());
	const Eigen::Index n  = static_cast<Eigen::Index>(r.NumNaturalFunctions());
	const Eigen::Index nv = static_cast<Eigen::Index>(r.NumVariables());

	Vec<complex_dbl> p(nv);
	for (Eigen::Index k = 0; k < nv; ++k) p(k) = complex_dbl(0.4 + 0.1 * static_cast<double>(k), 0.2 - 0.05 * static_cast<double>(k));

	// poison the output buffers with a huge value; the block must overwrite every owned entry.
	Vec<complex_dbl> seg(n);   seg.setConstant(complex_dbl(1e300, 1e300));
	Mat<complex_dbl> jac(n, nv); jac.setConstant(complex_dbl(1e300, 1e300));

	block.EvalInPlace<complex_dbl>(seg, p, complex_dbl(0));
	block.JacobianInPlace<complex_dbl>(jac, p, complex_dbl(0));

	// reference: the function-tree expansion of just the block's rows
	System twin = r.ExpandToFunctionTree();
	Vec<complex_dbl> seg_ref = twin.Eval(p);
	Mat<complex_dbl> jac_ref = twin.Jacobian(p);

	for (Eigen::Index i = 0; i < n; ++i)
		BOOST_CHECK(std::abs(seg(i) - seg_ref(i)) < 1e-10);
	for (Eigen::Index i = 0; i < n; ++i)
		for (Eigen::Index j = 0; j < nv; ++j)
			BOOST_CHECK(std::abs(jac(i, j) - jac_ref(i, j)) < 1e-10);
}


// ---- metadata ------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metadata)
{
	DefaultPrecision(30);
	System r = OverdeterminedSingleGroup().Randomize();
	auto const& block = std::get<Randomization>(r.Blocks().front());
	BOOST_CHECK_EQUAL(block.NumFunctions(), 2u);
	BOOST_CHECK(!block.DependsOnPathVariable());     // autonomous operand
	BOOST_CHECK(!block.HasConstantJacobian());

	// time-derivative is identically zero
	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	Vec<complex_dbl> dt(2); dt.setConstant(complex_dbl(7));
	block.TimeDerivInPlace<complex_dbl>(dt, p, complex_dbl(0));
	BOOST_CHECK(std::abs(dt(0)) < 1e-14 && std::abs(dt(1)) < 1e-14);
}

BOOST_AUTO_TEST_SUITE_END()
