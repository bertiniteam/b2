//This file is part of Bertini 2.
//
//moving_homotopy_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//moving_homotopy_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file moving_homotopy_test.cpp

\brief Unit tests for MakeMovingHomotopy: move only the moving rows, keep the fixed system as its
own blocks (evaluated once, contributing zero to dH/dt).
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/block.hpp"

BOOST_AUTO_TEST_SUITE(moving_homotopy_suite)

using namespace bertini;
using bertini::node::Variable;

// gamma fixed off the real axis for reproducibility.
static std::shared_ptr<node::Node> Gamma()
{
	return node::Complex::Make(mpfr_complex("0.6", "0.8"));
}

// A LinearFormsBlock-backed single-form system  c.[x;1]  over variables {x,y}.
static System LinearFormSystem(std::shared_ptr<node::Variable> x, std::shared_ptr<node::Variable> y,
                               mpfr_complex cx, mpfr_complex cy, mpfr_complex c1)
{
	System s;
	s.AddVariableGroup(VariableGroup{x, y});
	Mat<mpfr_complex> M(1, 3);
	M << cx, cy, c1;
	s.AddBlock(blocks::LinearFormsBlock(2, M));
	return s;
}

// H, its expansion twin, agree on values AND Jacobian at several (point, t), in dbl and mpfr.
static void AgreesWithExpansionAtTimes(System const& H)
{
	System twin = H.ExpandToFunctionTree();
	BOOST_REQUIRE_EQUAL(H.NumNaturalFunctions(), twin.NumNaturalFunctions());
	// (no Degrees() comparison: the blend reports space-only degree, while the expanded node tree
	//  counts the path variable t too, so they legitimately differ for a homotopy.)

	const Eigen::Index nv = static_cast<Eigen::Index>(H.NumVariables());
	std::vector<dbl> times{dbl(1.0), dbl(0.0), dbl(0.37, -0.21)};

	Vec<dbl> p(nv);
	for (Eigen::Index k = 0; k < nv; ++k) p(k) = dbl(0.3 + 0.13 * static_cast<double>(k), 0.4 - 0.07 * static_cast<double>(k));

	for (auto t : times)
	{
		Vec<dbl> a = H.Eval(p, t),  b = twin.Eval(p, t);
		for (Eigen::Index i = 0; i < a.size(); ++i)
			BOOST_CHECK(std::abs(a(i) - b(i)) < 1e-10);

		Mat<dbl> ja = H.Jacobian(p, t), jb = twin.Jacobian(p, t);
		for (Eigen::Index i = 0; i < ja.rows(); ++i)
			for (Eigen::Index j = 0; j < ja.cols(); ++j)
				BOOST_CHECK(std::abs(ja(i, j) - jb(i, j)) < 1e-10);

		Vec<dbl> da = H.TimeDerivative(p, t), db = twin.TimeDerivative(p, t);
		for (Eigen::Index i = 0; i < da.size(); ++i)
			BOOST_CHECK(std::abs(da(i) - db(i)) < 1e-10);
	}

	// one mpfr cross-check
	DefaultPrecision(40);
	Vec<mpfr_complex> pm(nv);
	for (Eigen::Index k = 0; k < nv; ++k) pm(k) = mpfr_complex(p(k).real(), p(k).imag());
	mpfr_complex tm("0.37", "-0.21");
	Vec<mpfr_complex> am = H.Eval(pm, tm), bm = twin.Eval(pm, tm);
	for (Eigen::Index i = 0; i < am.size(); ++i)
		BOOST_CHECK(abs(am(i) - bm(i)) < mpfr_float("1e-30"));
	Mat<mpfr_complex> jam = H.Jacobian(pm, tm), jbm = twin.Jacobian(pm, tm);
	for (Eigen::Index i = 0; i < jam.rows(); ++i)
		for (Eigen::Index j = 0; j < jam.cols(); ++j)
			BOOST_CHECK(abs(jam(i, j) - jbm(i, j)) < mpfr_float("1e-30"));
}


// fixed = unit circle (PolynomialBlock); moving slice y (start) -> y-x (end), both polynomial.
static System CircleMovingSlice()
{
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y});
	fixed.AddFunction(x*x + y*y - node::Integer::Make(1));
	System start_moving; start_moving.AddVariableGroup(VariableGroup{x, y}); start_moving.AddFunction(y);
	System end_moving;   end_moving.AddVariableGroup(VariableGroup{x, y});   end_moving.AddFunction(y - x);
	return MakeMovingHomotopy(fixed, start_moving, end_moving, "t", Gamma());
}


BOOST_AUTO_TEST_CASE(shape_and_path_variable)
{
	DefaultPrecision(30);
	System H = CircleMovingSlice();
	BOOST_CHECK_EQUAL(H.NumNaturalFunctions(), 2u);     // 1 fixed (circle) + 1 moving (slice)
	BOOST_CHECK(H.HavePathVariable());
}

BOOST_AUTO_TEST_CASE(endpoints_t0_and_t1)
{
	DefaultPrecision(40);
	System H = CircleMovingSlice();

	Vec<dbl> p(2); p << dbl(0.4, 0.2), dbl(-0.3, 0.5);
	const dbl circle = p(0)*p(0) + p(1)*p(1) - dbl(1);

	// t = 0: moving row = end_moving = y - x ; fixed row = circle.
	Vec<dbl> at0 = H.Eval(p, dbl(0));
	BOOST_CHECK(std::abs(at0(0) - circle) < 1e-12);
	BOOST_CHECK(std::abs(at0(1) - (p(1) - p(0))) < 1e-12);

	// t = 1: moving row = gamma * start_moving = gamma * y ; fixed row still circle.
	const dbl gamma(0.6, 0.8);
	Vec<dbl> at1 = H.Eval(p, dbl(1));
	BOOST_CHECK(std::abs(at1(0) - circle) < 1e-12);
	BOOST_CHECK(std::abs(at1(1) - gamma * p(1)) < 1e-12);
}

BOOST_AUTO_TEST_CASE(fixed_row_is_left_out_of_time_derivative)
{
	DefaultPrecision(40);
	System H = CircleMovingSlice();

	Vec<dbl> p(2); p << dbl(0.4, 0.2), dbl(-0.3, 0.5);
	// dH/dt = [ 0 (circle is t-independent) ; -end + gamma*start = -(y-x) + gamma*y ]
	Vec<dbl> dt = H.TimeDerivative(p, dbl(0.5));
	BOOST_CHECK(std::abs(dt(0)) < 1e-14);                                  // fixed row: exactly out
	const dbl gamma(0.6, 0.8);
	BOOST_CHECK(std::abs(dt(1) - (-(p(1) - p(0)) + gamma * p(1))) < 1e-12); // moving row
}

BOOST_AUTO_TEST_CASE(matches_expansion_polynomial_moving)
{
	DefaultPrecision(40);
	AgreesWithExpansionAtTimes(CircleMovingSlice());
}

BOOST_AUTO_TEST_CASE(static_slice_and_moving_slice_both_fixed)
{
	// fixed = circle (poly) + a STATIC linear-forms slice  x ; moving slice y -> y - x.
	DefaultPrecision(40);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y});
	fixed.AddFunction(x*x + y*y - node::Integer::Make(1));               // row 0: circle (PolynomialBlock)
	Mat<mpfr_complex> M(1, 3); M << mpfr_complex(1), mpfr_complex(0), mpfr_complex(0);
	fixed.AddBlock(blocks::LinearFormsBlock(2, M));                      // row 1: static slice x=0 (LinearFormsBlock)

	System start_moving = LinearFormSystem(x, y, mpfr_complex(0), mpfr_complex(1), mpfr_complex(0));   // y
	System end_moving;   end_moving.AddVariableGroup(VariableGroup{x, y}); end_moving.AddFunction(y - x);

	System H = MakeMovingHomotopy(fixed, start_moving, end_moving, "t", Gamma());
	BOOST_CHECK_EQUAL(H.NumNaturalFunctions(), 3u);

	// both fixed rows (circle, static slice) are out of dH/dt; only the moving row survives.
	Vec<dbl> p(2); p << dbl(0.4, 0.2), dbl(-0.3, 0.5);
	Vec<dbl> dt = H.TimeDerivative(p, dbl(0.5));
	BOOST_CHECK(std::abs(dt(0)) < 1e-14);
	BOOST_CHECK(std::abs(dt(1)) < 1e-14);
	BOOST_CHECK(std::abs(dt(2)) > 1e-3);

	AgreesWithExpansionAtTimes(H);
}

BOOST_AUTO_TEST_CASE(deform_products_of_linears_into_polynomial)
{
	// regeneration step: deform (x-1)(x+1) (a products-of-linears block) into the circle, with a
	// static slice y (= the fixed row).
	DefaultPrecision(40);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed = LinearFormSystem(x, y, mpfr_complex(0), mpfr_complex(1), mpfr_complex(0));  // y (static)

	System start_moving; start_moving.AddVariableGroup(VariableGroup{x, y});
	Mat<mpfr_complex> f(2, 3);
	f << mpfr_complex(1), mpfr_complex(0), mpfr_complex(-1),     // (x - 1)
	     mpfr_complex(1), mpfr_complex(0), mpfr_complex(1);      // (x + 1)
	start_moving.AddBlock(blocks::ProductsOfLinearsBlock(2, std::vector<Mat<mpfr_complex>>{f}));

	System end_moving; end_moving.AddVariableGroup(VariableGroup{x, y});
	end_moving.AddFunction(x*x + y*y - node::Integer::Make(1));

	System H = MakeMovingHomotopy(fixed, start_moving, end_moving, "t", Gamma());
	BOOST_CHECK_EQUAL(H.NumNaturalFunctions(), 2u);

	// at t=1 the moving row is gamma*(x-1)(x+1) = gamma*(x^2-1); at t=0 it is the circle.
	Vec<dbl> p(2); p << dbl(2.0), dbl(3.0);
	const dbl gamma(0.6, 0.8);
	Vec<dbl> at1 = H.Eval(p, dbl(1));
	BOOST_CHECK(std::abs(at1(1) - gamma * (p(0)*p(0) - dbl(1))) < 1e-10);
	Vec<dbl> at0 = H.Eval(p, dbl(0));
	BOOST_CHECK(std::abs(at0(1) - (p(0)*p(0) + p(1)*p(1) - dbl(1))) < 1e-10);

	// fixed (static slice) row out of dH/dt
	Vec<dbl> dt = H.TimeDerivative(p, dbl(0.5));
	BOOST_CHECK(std::abs(dt(0)) < 1e-14);

	AgreesWithExpansionAtTimes(H);
}

BOOST_AUTO_TEST_CASE(rejects_mismatched_endpoints)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y}); fixed.AddFunction(x*x + y*y - node::Integer::Make(1));
	System sm; sm.AddVariableGroup(VariableGroup{x, y}); sm.AddFunction(y);
	System em; em.AddVariableGroup(VariableGroup{x, y}); em.AddFunction(y - x); em.AddFunction(x);  // 2 != 1
	BOOST_CHECK_THROW(MakeMovingHomotopy(fixed, sm, em, "t", Gamma()), std::runtime_error);
}

// Regression (b2 issue #258): concatenating the fixed system INTO the moving rows (instead of
// passing only the rows that move) duplicates the fixed equations.  The count/structure checks pass
// (both moving systems have 2 functions over the same variables), so the duplication is caught by
// the structural top-level-function comparison instead.
BOOST_AUTO_TEST_CASE(rejects_fixed_function_duplicated_in_moving_rows)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	auto circle = x*x + y*y - node::Integer::Make(1);
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y}); fixed.AddFunction(circle);

	// the mistake: moving rows = fixed (circle) AND the slice, in both endpoints.
	System sm; sm.AddVariableGroup(VariableGroup{x, y}); sm.AddFunction(circle); sm.AddFunction(y);
	System em; em.AddVariableGroup(VariableGroup{x, y}); em.AddFunction(circle); em.AddFunction(y - x);
	BOOST_CHECK_THROW(MakeMovingHomotopy(fixed, sm, em, "t", Gamma()), std::runtime_error);

	// even a structurally-equal but independently-built copy of the fixed function is caught
	// (comparison is by serialized form, not pointer identity).
	System sm2; sm2.AddVariableGroup(VariableGroup{x, y}); sm2.AddFunction(x*x + y*y - node::Integer::Make(1)); sm2.AddFunction(y);
	System em2; em2.AddVariableGroup(VariableGroup{x, y}); em2.AddFunction(x*x + y*y - node::Integer::Make(1)); em2.AddFunction(y - x);
	BOOST_CHECK_THROW(MakeMovingHomotopy(fixed, sm2, em2, "t", Gamma()), std::runtime_error);

	// the correct call (slice-only moving rows) still builds.
	System sm_ok; sm_ok.AddVariableGroup(VariableGroup{x, y}); sm_ok.AddFunction(y);
	System em_ok; em_ok.AddVariableGroup(VariableGroup{x, y}); em_ok.AddFunction(y - x);
	BOOST_CHECK_NO_THROW(MakeMovingHomotopy(fixed, sm_ok, em_ok, "t", Gamma()));
}

// A "moving" row that is identical at both endpoints does not move; it belongs in `fixed`.  The
// check is positional (the blend pairs start/end rows by index), so row 1 here trips it while the
// genuinely-moving row 0 does not.
BOOST_AUTO_TEST_CASE(rejects_non_moving_row_in_moving_block)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x"), y = Variable::Make("y");
	System fixed; fixed.AddVariableGroup(VariableGroup{x, y}); fixed.AddFunction(x*x + y*y - node::Integer::Make(1));

	System sm; sm.AddVariableGroup(VariableGroup{x, y}); sm.AddFunction(y);     sm.AddFunction(x);  // row 1: x
	System em; em.AddVariableGroup(VariableGroup{x, y}); em.AddFunction(y - x); em.AddFunction(x);  // row 1: x (static!)
	BOOST_CHECK_THROW(MakeMovingHomotopy(fixed, sm, em, "t", Gamma()), std::runtime_error);
}

// Regression: Clone(System) is now a Memory-isolating shallow copy --- it shares the
// immutable node DAG and the compiled SLP Program, but copies the per-thread evaluation Memory at
// every level, INCLUDING a BlendBlock's nested operand Systems (which are deep-copied via the
// System copy constructor: own Memory, shared DAG).  A clone of a moving homotopy must therefore
// reproduce the original exactly, and the two must evaluate independently.  (Thread-safety itself
// --- no shared mutable state --- is by construction here and is exercised by the threaded zero-dim
// solves; a sequential test cannot observe a data race because each Eval re-sets its inputs.)
BOOST_AUTO_TEST_CASE(clone_of_moving_homotopy_reproduces_and_is_independent)
{
	DefaultPrecision(30);
	System H = CircleMovingSlice();
	H.Differentiate();

	Vec<dbl> p1(2); p1 << dbl(0.3, 0.1), dbl(-0.2, 0.4);
	Vec<dbl> p2(2); p2 << dbl(1.5, -0.7), dbl(0.9, 0.2);
	const dbl t1(0.25, 0.0), t2(0.8, -0.1);

	const Vec<dbl> f1 = H.Eval(p1, t1);
	const Mat<dbl> j1 = H.Jacobian(p1, t1);

	System H_clone = Clone(H);

	// The clone (its deep-copied operands included) reproduces the original exactly.
	BOOST_CHECK(H_clone.Eval(p1, t1).isApprox(f1));
	BOOST_CHECK(H_clone.Jacobian(p1, t1).isApprox(j1));

	// Evaluating the clone elsewhere leaves the original's results unchanged (independent state).
	(void) H_clone.Eval(p2, t2);
	(void) H_clone.Jacobian(p2, t2);
	BOOST_CHECK(H.Eval(p1, t1).isApprox(f1));
	BOOST_CHECK(H.Jacobian(p1, t1).isApprox(j1));
}

BOOST_AUTO_TEST_SUITE_END()
