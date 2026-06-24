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
}


// f0 = 2x + 3y + 1,  f1 = x - y + 4   (augmented rows: [coef_x, coef_y, const])
// at (x,y) = (1,1):  f0 = 6,  f1 = 4
static Slice MakeKnownSlice()
{
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<mpfr_complex> M(2,3);
	M << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1),
	     mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);

	return Slice::FromCoefficients(vars, M);
}


BOOST_AUTO_TEST_CASE(from_coefficients_eval)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	BOOST_CHECK_EQUAL(s.Dimension(),2);
	BOOST_CHECK_EQUAL(s.NumVariables(),2);

	Vec<dbl> p(2); p << dbl(1), dbl(1);
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

	Vec<dbl> p(2); p << dbl(1), dbl(1);
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
	BOOST_CHECK((h.Coefficients() - s.Coefficients().topRows(2)).norm() < mpfr_float("1e-30"));

	auto t = s.Tail(1);
	BOOST_CHECK_EQUAL(t.Dimension(),1);
	BOOST_CHECK((t.Coefficients() - s.Coefficients().bottomRows(1)).norm() < mpfr_float("1e-30"));

	auto r = s.Rows(std::vector<unsigned>{0,3});
	BOOST_CHECK_EQUAL(r.Dimension(),2);
	BOOST_CHECK((r.Coefficients().row(0) - s.Coefficients().row(0)).norm() < mpfr_float("1e-30"));
	BOOST_CHECK((r.Coefficients().row(1) - s.Coefficients().row(3)).norm() < mpfr_float("1e-30"));

	BOOST_CHECK_THROW(s.Head(99), std::runtime_error);
	BOOST_CHECK_THROW(s.Rows(std::vector<unsigned>{99}), std::runtime_error);
}


BOOST_AUTO_TEST_CASE(homogeneous_has_zero_constant_column)
{
	Var x = Variable::Make("x"), y = Variable::Make("y"), z = Variable::Make("z");
	VariableGroup vars{x,y,z};

	auto s = Slice::RandomComplex(vars, 2, /*homogeneous=*/true);

	BOOST_CHECK(s.IsHomogeneous());
	BOOST_CHECK(s.Coefficients().rightCols(1).norm() < mpfr_float("1e-30"));
}


BOOST_AUTO_TEST_CASE(precision_roundtrip)
{
	DefaultPrecision(30);
	auto s = MakeKnownSlice();

	s.Precision(50);
	BOOST_CHECK_EQUAL(s.Precision(),50u);

	Vec<mpfr_complex> p(2); p << mpfr_complex(1), mpfr_complex(1);
	auto v = s.Eval(p);
	BOOST_CHECK(abs(v(0) - mpfr_complex(6)) < mpfr_float("1e-25"));
	BOOST_CHECK(abs(v(1) - mpfr_complex(4)) < mpfr_float("1e-25"));
}


// A slice added to a System evaluates identically through the System.
BOOST_AUTO_TEST_CASE(add_to_system_agrees_with_slice_eval)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<mpfr_complex> M(2,3);
	M << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1),
	     mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);
	auto s = Slice::FromCoefficients(vars, M);

	System sys;
	sys.AddVariableGroup(vars);
	s.AddTo(sys);

	BOOST_CHECK_EQUAL(sys.NumNaturalFunctions(), 2u);

	Vec<dbl> p(2); p << dbl(1), dbl(1);
	auto from_system = sys.Eval(p);
	auto from_slice  = s.Eval(p);

	BOOST_CHECK_EQUAL(from_system.size(), from_slice.size());
	BOOST_CHECK_SMALL((from_system - from_slice).norm(), 1e-11);
}


BOOST_AUTO_TEST_CASE(concatenate_stacks_forms)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y");
	VariableGroup vars{x,y};

	Mat<mpfr_complex> A(1,3); A << mpfr_complex(2), mpfr_complex(3),  mpfr_complex(1);  // 2x + 3y + 1
	Mat<mpfr_complex> B(1,3); B << mpfr_complex(1), mpfr_complex(-1), mpfr_complex(4);  // x - y + 4
	auto sa = Slice::FromCoefficients(vars, A);
	auto sb = Slice::FromCoefficients(vars, B);

	auto s = sa.Concatenate(sb);
	BOOST_CHECK_EQUAL(s.Dimension(), 2);
	BOOST_CHECK_EQUAL(s.NumVariables(), 2);

	Vec<dbl> p(2); p << dbl(1), dbl(1);
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

	Vec<dbl> p(2); p << dbl(1), dbl(1);
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
	Vec<dbl> p(2); p << dbl(1), dbl(1);
	auto v = s2.Eval(p);
	BOOST_CHECK_CLOSE(v(0).real(), 6.0, 1e-11);
	BOOST_CHECK_CLOSE(v(1).real(), 4.0, 1e-11);
}


BOOST_AUTO_TEST_SUITE_END()
