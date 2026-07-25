//This file is part of Bertini 2.
//
//test/nag_algorithms/newton_refine.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/nag_algorithms/newton_refine.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/nag_algorithms/newton_refine.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// silviana amethyst, university of wisconsin eau claire

/**
\file test/nag_algorithms/newton_refine.cpp

Tests for the standalone NewtonRefine -- above all, refinement of SINGULAR points on
DEFLATED systems, which is the whole point: deflation restores quadratic convergence
exactly where the plain system cannot converge.  The specimens walk the isosingular
hierarchy: the double cone's isolated singularity (one deflation), a non-origin point
of the whitney umbrella's handle (a smooth point of the singular curve; one deflation
plus a pinning slice), and the umbrella's origin (the pinch point -- a singular point
OF the singular embedded curve, needing the SECOND deflation stage).

Overdetermined deflated systems are squared by randomization with hardcoded
small-prime matrices whose nonsingularity at the target point is verified by hand in
comments -- the same recipe the tracking layer uses in production.
*/

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "bertini2/nag_algorithms/newton_refine.hpp"
#include "bertini2/system/system.hpp"

using bertini::System;
using bertini::Vec;
using bertini::Mat;
using bertini::complex_mp;
using bertini::DefaultPrecision;
using bertini::node::Variable;
using bertini::algorithm::NewtonRefine;
using bertini::SuccessCode;

BOOST_AUTO_TEST_SUITE(standalone_newton_refine)


// a nonsingular root: full quadratic convergence, achieved accuracy at the ask
BOOST_AUTO_TEST_CASE(nonsingular_root_refines_quadratically)
{
	DefaultPrecision(60);
	auto x = Variable::Make("x");
	System S;
	S.AddUngroupedVariable(x);
	S.AddFunction(pow(x,2) - 2);

	Vec<complex_mp> start(1);
	start << complex_mp("1.4142");           // ~1.4e-5 from sqrt(2)

	auto r = NewtonRefine(S, start, 1e-40, 50);
	BOOST_CHECK(r.code == SuccessCode::Success);
	BOOST_CHECK(r.achieved <= 1e-40);
	BOOST_CHECK(r.iterations <= 10);         // quadratic: ~3 doublings needed
	using mpfr_float = bertini::real_mp;
	mpfr_float residual = abs(pow(r.point(0),2) - complex_mp(2));
	BOOST_CHECK(residual < mpfr_float("1e-38"));
}


// overdetermined systems are refused with instructions, not mangled
BOOST_AUTO_TEST_CASE(overdetermined_system_is_refused)
{
	DefaultPrecision(30);
	auto x = Variable::Make("x");
	System S;
	S.AddUngroupedVariable(x);
	S.AddFunction(pow(x,2));
	S.AddFunction(x - 1);

	Vec<complex_mp> start(1);
	start << complex_mp("0.5");
	BOOST_CHECK_THROW(NewtonRefine(S, start, 1e-20, 10), std::runtime_error);
}


// the disease, isolated: at a double root the plain system converges only
// linearly (steps halve), so a tight tolerance is out of reach in few iterations
BOOST_AUTO_TEST_CASE(plain_newton_stalls_at_a_double_root)
{
	DefaultPrecision(60);
	auto x = Variable::Make("x");
	System S;
	S.AddUngroupedVariable(x);
	S.AddFunction(pow(x,2));

	Vec<complex_mp> start(1);
	start << complex_mp("1e-6");
	auto r = NewtonRefine(S, start, 1e-40, 25);
	BOOST_CHECK(r.code == SuccessCode::FailedToConverge);
	BOOST_CHECK(r.achieved > 1e-40);         // ~1e-6/2^25 ~ 3e-14: nowhere near
}


// silviana specimen 1: the double cone x^2+y^2-z^2, singular at the origin.
// Deflation appends the gradient; the overdetermined [f; grad f] (4 fns, 3 vars)
// is squared by the hardcoded randomization
//     R = [ 1  2  3  5 ]
//         [ 7 11 13 17 ]
//         [19 23 29 31 ]
// At the origin grad f = 0 and J([f;fx;fy;fz]) has rows 0,(2,0,0),(0,2,0),(0,0,-2),
// so J(R.F)(0) has rows 2*(2,3,-5), 2*(11,13,-17), 2*(23,29,-31) with
// det(base) = -70 != 0: the deflated randomized system is REGULAR at the origin.
BOOST_AUTO_TEST_CASE(double_cone_singularity_refines_on_deflated_system)
{
	DefaultPrecision(60);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	auto f  = pow(x,2) + pow(y,2) - pow(z,2);
	auto fx = 2*x;
	auto fy = 2*y;
	auto fz = -2*z;

	System S;
	bertini::VariableGroup vg{x, y, z};
	S.AddVariableGroup(vg);
	S.AddFunction( 1*f +  2*fx +  3*fy +  5*fz);
	S.AddFunction( 7*f + 11*fx + 13*fy + 17*fz);
	S.AddFunction(19*f + 23*fx + 29*fy + 31*fz);

	Vec<complex_mp> start(3);
	start << complex_mp("1e-6"), complex_mp("-2e-6"), complex_mp("5e-7");

	auto r = NewtonRefine(S, start, 1e-45, 50);
	BOOST_CHECK(r.code == SuccessCode::Success);
	BOOST_CHECK(r.achieved <= 1e-45);
	using mpfr_float = bertini::real_mp;
	mpfr_float dist = max(abs(r.point(0)), max(abs(r.point(1)), abs(r.point(2))));
	BOOST_CHECK(dist < mpfr_float("1e-40"));  // landed ON the singularity
}


// silviana specimen 2: a NON-ORIGIN point of the whitney umbrella's handle.
// f = x^2 - y^2 z is singular along the whole handle x=y=0; the point (0,0,1) is a
// SMOOTH point of that singular curve.  One deflation (the gradient
// {2x, -2yz, -y^2}) plus the pinning slice z-1 gives 4 fns over 3 vars, squared by
//     R = [ 1  2  3  5 ]
//         [ 7 11 13 17 ]
//         [19 23 29 31 ]
// At (0,0,1) the stacked Jacobian has rows (2,0,0),(0,-2,0),(0,0,0),(0,0,1), so
// J(R.F) has rows (2,-4,5),(14,-22,17),(38,-46,31) with det = 312 != 0: REGULAR.
BOOST_AUTO_TEST_CASE(whitney_handle_point_refines_on_deflated_system)
{
	DefaultPrecision(60);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	auto g1 = 2*x;              // f_x
	auto g2 = -2*y*z;           // f_y
	auto g3 = -pow(y,2);        // f_z
	auto g4 = z - 1;            // the pinning slice through the target point

	System S;
	bertini::VariableGroup vg{x, y, z};
	S.AddVariableGroup(vg);
	S.AddFunction( 1*g1 +  2*g2 +  3*g3 +  5*g4);
	S.AddFunction( 7*g1 + 11*g2 + 13*g3 + 17*g4);
	S.AddFunction(19*g1 + 23*g2 + 29*g3 + 31*g4);

	Vec<complex_mp> start(3);
	start << complex_mp("1e-7"), complex_mp("-1e-7"),
	         complex_mp(1) + complex_mp("1e-7");

	auto r = NewtonRefine(S, start, 1e-45, 50);
	BOOST_CHECK(r.code == SuccessCode::Success);
	BOOST_CHECK(r.achieved <= 1e-45);
	using mpfr_float = bertini::real_mp;
	mpfr_float dist = max(abs(r.point(0)),
	                      max(abs(r.point(1)), abs(r.point(2) - complex_mp(1))));
	BOOST_CHECK(dist < mpfr_float("1e-40"));  // landed ON the handle point
}


// silviana specimen 3: the whitney umbrella's ORIGIN -- the pinch point, a singular
// point OF the singular embedded curve.  The first deflation F1 = {2x, -2yz, -y^2}
// is itself singular there (its Jacobian has rank 1 at 0), so the SECOND stage
// appends 2x2 minors of J(F1): det[[2,0],[0,-2z]] = -4z and the row-1/row-3
// minor -4y.  F2 = {2x, -2yz, -y^2, -4z, -4y} (5 fns, 3 vars), squared by
//     R = [ 1  2  3  5  7 ]
//         [11 13 17 19 23 ]
//         [29 31 37 41 43 ]
// At the origin the stacked Jacobian has rows (2,0,0),0,0,(0,0,-4),(0,-4,0), so
// J(R.F2)(0) has rows (2,-28,-20),(22,-92,-76),(58,-172,-164); factoring 2 from
// column 1, det(base) = -2304 != 0: the SECOND-stage system is REGULAR.
BOOST_AUTO_TEST_CASE(whitney_pinch_point_needs_and_gets_second_deflation)
{
	DefaultPrecision(60);
	auto x = Variable::Make("x");
	auto y = Variable::Make("y");
	auto z = Variable::Make("z");
	auto g1 = 2*x;
	auto g2 = -2*y*z;
	auto g3 = -pow(y,2);
	auto m1 = -4*z;             // second-stage minor
	auto m2 = -4*y;             // second-stage minor

	// control: the FIRST-stage deflation alone, square as-is, is still singular at
	// the pinch point -- Newton limps and cannot reach a tight tolerance
	{
		System S1;
		bertini::VariableGroup vg{x, y, z};
		S1.AddVariableGroup(vg);
		S1.AddFunction(g1);
		S1.AddFunction(g2);
		S1.AddFunction(g3);
		Vec<complex_mp> start(3);
		start << complex_mp("1e-7"), complex_mp("1e-7"), complex_mp("1e-7");
		auto r1 = NewtonRefine(S1, start, 1e-45, 25);
		BOOST_CHECK(r1.code != SuccessCode::Success);
	}

	// the second-stage deflated randomized system nails it
	System S2;
	bertini::VariableGroup vg{x, y, z};
	S2.AddVariableGroup(vg);
	S2.AddFunction( 1*g1 +  2*g2 +  3*g3 +  5*m1 +  7*m2);
	S2.AddFunction(11*g1 + 13*g2 + 17*g3 + 19*m1 + 23*m2);
	S2.AddFunction(29*g1 + 31*g2 + 37*g3 + 41*m1 + 43*m2);

	Vec<complex_mp> start(3);
	start << complex_mp("1e-7"), complex_mp("1e-7"), complex_mp("1e-7");

	auto r = NewtonRefine(S2, start, 1e-45, 50);
	BOOST_CHECK(r.code == SuccessCode::Success);
	BOOST_CHECK(r.achieved <= 1e-45);
	using mpfr_float = bertini::real_mp;
	mpfr_float dist = max(abs(r.point(0)), max(abs(r.point(1)), abs(r.point(2))));
	BOOST_CHECK(dist < mpfr_float("1e-40"));  // landed ON the pinch point
}


BOOST_AUTO_TEST_SUITE_END()
