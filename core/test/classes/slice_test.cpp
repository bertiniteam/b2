//This file is part of Bertini 2.
//
//slice_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//slice_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with slice_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

//  slice_test.cpp
//

/**
\file slice_test.cpp Unit testing for slicing
*/

#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"

BOOST_AUTO_TEST_SUITE(linear_slicing)

using namespace bertini;
using Var = std::shared_ptr<bertini::node::Variable>;
using Variable = bertini::node::Variable;

using bertini::DefaultPrecision;

BOOST_AUTO_TEST_CASE(slice_basic_complex)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	auto s = Slice::RandomComplex(vars,2);

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);
	// the augmented coefficient matrix is (#forms) x (#vars + 1)
	BOOST_CHECK_EQUAL(s.Coefficients().rows(),2);
	BOOST_CHECK_EQUAL(s.Coefficients().cols(),4);
}


BOOST_AUTO_TEST_CASE(slice_basic_crazy_overslice)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	auto s = Slice::RandomComplex(vars,6);

	BOOST_CHECK_EQUAL(s.Dimension(),6);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);
}


BOOST_AUTO_TEST_CASE(slice_basic_real)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	VariableGroup vars{x,y,z};

	auto s = Slice::RandomReal(vars,2);

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),3);

	// issue #294: a real slice's coefficients must actually be real.  Previously the orthogonal
	// (default) path fell through to the COMPLEX conjugate-orthonormal matrix, so every "real"
	// slice came out complex.  Every entry of the augmented coefficient matrix has zero imaginary part.
	auto const& C = s.Coefficients();
	for (int ii = 0; ii < C.rows(); ++ii)
		for (int jj = 0; jj < C.cols(); ++jj)
			BOOST_CHECK_EQUAL(C(ii,jj).imag(), real_mp(0));
}


// f0 = 2x + 3y + 1,  f1 = x - y + 4   (augmented rows: [coef_x, coef_y, const])
// at (x,y) = (1,1):  f0 = 6,  f1 = 4
static Slice MakeKnownSlice()
{
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<complex_mp> M(2,3);
	M << complex_mp(2), complex_mp(3),  complex_mp(1),
	     complex_mp(1), complex_mp(-1), complex_mp(4);

	return Slice::FromCoefficients(vars, M);
}


BOOST_AUTO_TEST_CASE(from_coefficients_eval)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),2);

	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto v = s.Eval(p);

	BOOST_CHECK_EQUAL(v.size(),2);
	BOOST_CHECK_CLOSE(v(0).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(v(1).real(), 4.0, 1e-11);
	BOOST_CHECK_SMALL(v(0).imag(), 1e-11);
	BOOST_CHECK_SMALL(v(1).imag(), 1e-11);
}


BOOST_AUTO_TEST_CASE(jacobian_is_variable_coefficient_block)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto J = s.Jacobian(p);

	BOOST_CHECK_EQUAL(J.rows(),2);
	BOOST_CHECK_EQUAL(J.cols(),2);   // the constant column is dropped
	BOOST_CHECK_CLOSE(J(0,0).real(),  2.0, 1e-11);
	BOOST_CHECK_CLOSE(J(0,1).real(),  3.0, 1e-11);
	BOOST_CHECK_CLOSE(J(1,0).real(),  1.0, 1e-11);
	BOOST_CHECK_CLOSE(J(1,1).real(), -1.0, 1e-11);
}


BOOST_AUTO_TEST_CASE(head_tail_rows_subsetting)
{
	DefaultPrecision(30);
	Var w = Variable::Make("w"), x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{w,x,y};

	auto s = Slice::RandomComplex(vars,4);

	auto h = s.Head(2);
	BOOST_CHECK_EQUAL(h.Dimension(),2);
	BOOST_CHECK_EQUAL(h.NumVariables(),3);
	// the head's forms are exactly the slice's first two forms.
	BOOST_CHECK((h.Coefficients() - s.Coefficients().topRows(2)).norm() < real_mp("1e-30"));

	auto t = s.Tail(1);
	BOOST_CHECK_EQUAL(t.Dimension(),1);
	BOOST_CHECK((t.Coefficients() - s.Coefficients().bottomRows(1)).norm() < real_mp("1e-30"));

	auto r = s.Rows(std::vector<unsigned>{0,3});
	BOOST_CHECK_EQUAL(r.Dimension(),2);
	BOOST_CHECK((r.Coefficients().row(0) - s.Coefficients().row(0)).norm() < real_mp("1e-30"));
	BOOST_CHECK((r.Coefficients().row(1) - s.Coefficients().row(3)).norm() < real_mp("1e-30"));

	BOOST_CHECK_THROW(s.Head(99), std::runtime_error);
	BOOST_CHECK_THROW(s.Rows(std::vector<unsigned>{99}), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(homogeneous_has_zero_constant_column)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	auto s = Slice::RandomComplex(vars, 2, /*homogeneous=*/true);

	BOOST_CHECK(s.IsHomogeneous());
	BOOST_CHECK(s.Coefficients().rightCols(1).norm() < real_mp("1e-30"));
}


BOOST_AUTO_TEST_CASE(precision_roundtrip)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	s.Precision(50);
	BOOST_CHECK_EQUAL(s.Precision(),50u);

	Vec<complex_mp> p(2); p << complex_mp(1), complex_mp(1);
	auto v = s.Eval(p);
	BOOST_CHECK(abs(v(0) - complex_mp(6)) < real_mp("1e-25"));
	BOOST_CHECK(abs(v(1) - complex_mp(4)) < real_mp("1e-25"));
}


// A slice added to a System evaluates identically through the System.
BOOST_AUTO_TEST_CASE(add_to_system_agrees_with_slice_eval)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<complex_mp> M(2,3);
	M << complex_mp(2), complex_mp(3),  complex_mp(1),
	     complex_mp(1), complex_mp(-1), complex_mp(4);
	auto s = Slice::FromCoefficients(vars, M);

	System sys;
	sys.AddVariableGroup(vars);
	s.AddTo(sys);

	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2u);

	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto from_system = sys.Eval(p);
	auto from_slice  = s.Eval(p);

	BOOST_CHECK_EQUAL(from_system.size(), from_slice.size());
	BOOST_CHECK_SMALL((from_system - from_slice).norm(), 1e-11);
}


BOOST_AUTO_TEST_CASE(add_to_homogenized_system_folds_constant_onto_hom_var)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	// an affine system, homogenized: it now has a homogenizing variable.
	System sys;
	sys.AddVariableGroup(vars);
	sys.AddFunction(x*x + y*y - 1);
	sys.Homogenize();
	BOOST_CHECK_EQUAL(sys.NumHomVariables(), 1u);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 3u);

	// an affine slice (built on the two natural variables) added to the homogenized system: AddTo
	// folds its constant onto the homogenizing variable, so the added form is homogeneous degree 1.
	Mat<complex_mp> M(1,3); M << complex_mp(2), complex_mp(3), complex_mp(1); // 2x + 3y + 1
	auto s = Slice::FromCoefficients(vars, M);
	s.AddTo(sys);

	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2u);
	BOOST_CHECK(sys.IsHomogeneous());   // the whole system, slice form included, is homogeneous

	// the slice form is now a*x + b*y + c*h with c the old constant (1): vanishes at the origin.
	Vec<complex_dbl> origin(3); origin << complex_dbl(0), complex_dbl(0), complex_dbl(0);
	auto v = sys.Eval(origin);
	BOOST_CHECK_EQUAL(v.size(), 2);
	BOOST_CHECK_SMALL(std::abs(v(1)), 1e-11);   // the homogeneous slice form is 0 at the origin
}


BOOST_AUTO_TEST_CASE(add_to_system_rejects_variable_count_mismatch)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");

	// a slice on three variables cannot be added to a system on two.
	auto s = Slice::RandomComplex(VariableGroup{x,y,z}, 1);

	System sys;
	sys.AddVariableGroup(VariableGroup{x,y});

	BOOST_CHECK_THROW(s.AddTo(sys), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(concatenate_stacks_forms)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<complex_mp> A(1,3); A << complex_mp(2), complex_mp(3),  complex_mp(1);  // 2x + 3y + 1
	Mat<complex_mp> B(1,3); B << complex_mp(1), complex_mp(-1), complex_mp(4);  // x - y + 4
	auto sa = Slice::FromCoefficients(vars, A);
	auto sb = Slice::FromCoefficients(vars, B);

	auto s = sa.Concatenate(sb);
	BOOST_CHECK_EQUAL(s.Dimension(), 2);
	BOOST_CHECK_EQUAL(s.NumVariables(), 2);

	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto v = s.Eval(p);
	BOOST_CHECK_CLOSE(v(0).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(v(1).real(), 4.0, 1e-11);

	// slices on different numbers of variables cannot be concatenated.
	Var z = Variable::Make("z");
	auto sc = Slice::RandomComplex(VariableGroup{x,y,z}, 1);
	BOOST_CHECK_THROW(sa.Concatenate(sc), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(as_system_evaluates_like_slice)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	auto sys = s.AsSystem();
	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2);
	BOOST_CHECK_EQUAL(sys.NumVariables(), 2);

	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	BOOST_CHECK_SMALL((sys.Eval(p) - s.Eval(p)).norm(), 1e-11);
}


BOOST_AUTO_TEST_CASE(serialization_roundtrip)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	std::stringstream ss;
	{ boost::archive::text_oarchive oa(ss); oa << s; }

	Slice s2;
	{ boost::archive::text_iarchive ia(ss); ia >> s2; }

	BOOST_CHECK_EQUAL(s2.Dimension(), 2);
	BOOST_CHECK_EQUAL(s2.NumVariables(), 2);
	BOOST_CHECK(!s2.IsHomogeneous());

	// the deserialized slice evaluates exactly like the original.
	Vec<complex_dbl> p(2); p << complex_dbl(1), complex_dbl(1);
	auto v = s2.Eval(p);
	BOOST_CHECK_CLOSE(v(0).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(v(1).real(), 4.0, 1e-11);
}


// ---- through a point (affine) --------------------------------------------------------

// The exact primitive: given a bare coefficient block A and a point p, the augmented matrix
// is [A | -A*p], so every form vanishes at p.  A = [[2,3],[1,-1]], p = (2,-3) gives constants
// -(2*2+3*-3)=5 and -(1*2-1*-3)=-5.
BOOST_AUTO_TEST_CASE(through_point_primitive_vanishes)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<complex_mp> A(2,2);
	A << complex_mp(2), complex_mp(3),
	     complex_mp(1), complex_mp(-1);
	Vec<complex_mp> p(2); p << complex_mp(2), complex_mp(-3);

	auto s = Slice::ThroughPoint(vars, A, p);

	BOOST_CHECK_EQUAL(s.Dimension(), 2);
	BOOST_CHECK(!s.IsHomogeneous());
	// the left (variable-coefficient) block is exactly A
	BOOST_CHECK(abs(s.Coefficients()(0,0) - complex_mp(2))  < real_mp("1e-25"));
	BOOST_CHECK(abs(s.Coefficients()(1,1) - complex_mp(-1)) < real_mp("1e-25"));
	// the assembled constant column is -A*p = (5, -5)
	BOOST_CHECK(abs(s.Coefficients()(0,2) - complex_mp(5))  < real_mp("1e-25"));
	BOOST_CHECK(abs(s.Coefficients()(1,2) - complex_mp(-5)) < real_mp("1e-25"));

	// every form vanishes at p
	Vec<complex_dbl> p_dbl(2); p_dbl << complex_dbl(2), complex_dbl(-3);
	BOOST_CHECK_SMALL(s.Eval(p_dbl).norm(), 1e-11);
}


BOOST_AUTO_TEST_CASE(through_point_rejects_wrong_length)
{
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<complex_mp> A(2,2);
	A << complex_mp(2), complex_mp(3),
	     complex_mp(1), complex_mp(-1);

	// point of the wrong length
	Vec<complex_mp> p3(3); p3 << complex_mp(2), complex_mp(-3), complex_mp(7);
	BOOST_CHECK_THROW(Slice::ThroughPoint(vars, A, p3), std::runtime_error);

	// coefficient block with the wrong number of columns (not num_variables)
	Mat<complex_mp> Awide(1,3); Awide << complex_mp(2), complex_mp(3), complex_mp(4);
	Vec<complex_mp> p2(2); p2 << complex_mp(2), complex_mp(-3);
	BOOST_CHECK_THROW(Slice::ThroughPoint(vars, Awide, p2), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(random_complex_through_point_vanishes)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	Vec<complex_mp> p(3); p << complex_mp(2), complex_mp(-3), complex_mp(5);

	auto s = Slice::RandomComplex(vars, 2, /*homogeneous=*/false, /*orthogonal=*/true, &p);

	BOOST_CHECK_EQUAL(s.Dimension(), 2);
	BOOST_CHECK(!s.IsHomogeneous());
	// the coefficient block is non-degenerate (not all zero)
	BOOST_CHECK(s.Coefficients().leftCols(3).norm() > real_mp("1e-3"));

	Vec<complex_dbl> p_dbl(3); p_dbl << complex_dbl(2), complex_dbl(-3), complex_dbl(5);
	BOOST_CHECK_SMALL(s.Eval(p_dbl).norm(), 1e-10);
}


BOOST_AUTO_TEST_CASE(random_real_through_point_vanishes)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	Vec<complex_mp> p(3); p << complex_mp(2), complex_mp(-3), complex_mp(5);

	auto s = Slice::RandomReal(vars, 2, /*homogeneous=*/false, /*orthogonal=*/true, &p);

	BOOST_CHECK(!s.IsHomogeneous());

	Vec<complex_dbl> p_dbl(3); p_dbl << complex_dbl(2), complex_dbl(-3), complex_dbl(5);
	BOOST_CHECK_SMALL(s.Eval(p_dbl).norm(), 1e-10);
}


// A real slice through a (real) point keeps a real variable-coefficient block; only the
// constant column carries -A*p (which would be complex only for a complex point).
BOOST_AUTO_TEST_CASE(real_through_point_real_coefficients)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	Vec<complex_mp> p(3); p << complex_mp(2), complex_mp(-3), complex_mp(5);

	auto s = Slice::RandomReal(vars, 2, /*homogeneous=*/false, /*orthogonal=*/true, &p);

	auto const& C = s.Coefficients();
	for (int ii = 0; ii < C.rows(); ++ii)
		for (int jj = 0; jj < C.cols() - 1; ++jj)   // the variable block (all but the constant column)
			BOOST_CHECK_EQUAL(C(ii,jj).imag(), real_mp(0));
}


// Temporary guard (task 6 replaces this with the homogeneous through-point construction):
// a slice cannot yet be both homogeneous and through a point.
BOOST_AUTO_TEST_CASE(through_point_rejects_homogeneous)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	Vec<complex_mp> p(3); p << complex_mp(2), complex_mp(-3), complex_mp(5);
	BOOST_CHECK_THROW(Slice::RandomComplex(vars, 2, /*homogeneous=*/true, /*orthogonal=*/true, &p),
	                  std::runtime_error);
}


BOOST_AUTO_TEST_SUITE_END()
