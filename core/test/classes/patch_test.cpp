//This file is part of Bertini 2.
//
//patch_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//patch_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with patch_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

//  patch_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.


/**
\file patch_test.cpp Unit testing for the bertini::Patch class.
*/


#include <boost/test/unit_test.hpp>
#include "bertini2/system/patch.hpp"




#include "externs.hpp"




BOOST_AUTO_TEST_SUITE(patch_class)


using bertini::DefaultPrecision;

template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;
using Patch = bertini::Patch;

using complex_dbl = bertini::complex_dbl;
using mpfr = bertini::complex_mp;

using real_mp = bertini::real_mp;


BOOST_AUTO_TEST_CASE(patch_create)
{
	Patch p;
}


BOOST_AUTO_TEST_CASE(patch_create_one_variable_group)
{
	std::vector<unsigned> s{2};

	Patch p(s);
}


BOOST_AUTO_TEST_CASE(patch_create_two_variable_groups)
{
	std::vector<unsigned> s{2,3};

	Patch p(s);
}


BOOST_AUTO_TEST_CASE(patch_eval_two_variable_groups_prec16)
{
	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<complex_dbl> v(5);
	v << complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1);

	p.Precision(16);
	
	auto f = p.Eval(v);

	BOOST_CHECK_EQUAL(f.size(),2);
}


BOOST_AUTO_TEST_CASE(patch_jacobian_two_variable_groups_prec16)
{
	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<complex_dbl> v(5);
	v << complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1);

	p.Precision(16);

	auto J = p.Jacobian(v);
	BOOST_CHECK_EQUAL(J.rows(),2);
	BOOST_CHECK_EQUAL(J.cols(),5);

	BOOST_CHECK_EQUAL(J(0,2),complex_dbl(0));
	BOOST_CHECK_EQUAL(J(0,3),complex_dbl(0));
	BOOST_CHECK_EQUAL(J(0,4),complex_dbl(0));

	BOOST_CHECK_EQUAL(J(1,0),complex_dbl(0));
	BOOST_CHECK_EQUAL(J(1,1),complex_dbl(0));
}


// A patch row is sparse (one coefficient block per variable group, zero elsewhere), so
// JacobianInPlace must FULLY define the rows it owns -- including zeroing the off-coefficient
// entries -- rather than relying on the caller to hand it a zeroed matrix.  The block-composed
// System::Jacobian path allocates J uninitialized and assigns only the function (block) rows, so a
// patch that left its off-coefficient entries untouched read uninitialized memory there: benign
// (zeroed pages) on Linux/macOS, but garbage on Windows, where a degree-2 homotopy's Jacobian
// "evaluated" to ~1e252 and wrecked the AMP condition-number estimate.  Hand it a fully-poisoned
// buffer and require every owned entry to be correct.
BOOST_AUTO_TEST_CASE(patch_jacobian_fully_defines_its_rows_into_a_dirty_buffer)
{
	std::vector<unsigned> s{2,3};

	Patch p(s);
	p.Precision(16);

	Vec<complex_dbl> v(5);
	v << complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1);

	Mat<complex_dbl> J = Mat<complex_dbl>::Constant(2, 5, complex_dbl(1e300)); // poison every entry

	p.JacobianInPlace(J, v);

	// off-coefficient entries of each patch row must be overwritten with zero, not left poisoned
	BOOST_CHECK_EQUAL(J(0,2), complex_dbl(0));
	BOOST_CHECK_EQUAL(J(0,3), complex_dbl(0));
	BOOST_CHECK_EQUAL(J(0,4), complex_dbl(0));
	BOOST_CHECK_EQUAL(J(1,0), complex_dbl(0));
	BOOST_CHECK_EQUAL(J(1,1), complex_dbl(0));

	// the coefficient entries are still written (no longer the poison value)
	BOOST_CHECK_NE(J(0,0), complex_dbl(1e300));
	BOOST_CHECK_NE(J(1,2), complex_dbl(1e300));
}



BOOST_AUTO_TEST_CASE(patch_eval_two_variable_groups_prec30)
{
	DefaultPrecision(30);
	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<mpfr> v(5);
	v << mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1);

	p.Precision(30);
	auto f = p.Eval(v);
}


BOOST_AUTO_TEST_CASE(patch_jacobian_two_variable_groups_prec30)
{
	DefaultPrecision(30);

	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<mpfr> v(5);
	v << mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1);

	p.Precision(30);

	auto J = p.Jacobian(v);

	BOOST_CHECK_EQUAL(J.rows(),2);
	BOOST_CHECK_EQUAL(J.cols(),5);

	BOOST_CHECK_EQUAL(J(0,2),mpfr(0));
	BOOST_CHECK_EQUAL(J(0,3),mpfr(0));
	BOOST_CHECK_EQUAL(J(0,4),mpfr(0));

	BOOST_CHECK_EQUAL(J(1,0),mpfr(0));
	BOOST_CHECK_EQUAL(J(1,1),mpfr(0));
}





BOOST_AUTO_TEST_CASE(patch_rescale_and_evaluate_prec16)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::vector<unsigned> s{2,3};

	Patch p(s);
	p.Precision(16);

	Vec<complex_dbl> v(5);
	v << complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1),  complex_dbl(1);

	auto v_rescaled = p.RescalePoint(v);

	auto f = p.Eval(v_rescaled);

	BOOST_CHECK_EQUAL(f.size(),2);
	for (int ii = 0; ii < 2; ++ii)
		BOOST_CHECK(abs(f(ii)) < threshold_clearance_d);
	
}




BOOST_AUTO_TEST_CASE(patch_rescale_and_evaluate_prec_default_mpfr)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<mpfr> v(5);
	v << mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1),  mpfr(1);

	auto v_rescaled = p.RescalePoint(v);

	auto f = p.Eval(v_rescaled);

	BOOST_CHECK_EQUAL(f.size(),2);
	for (int ii = 0; ii < 2; ++ii)
		BOOST_CHECK(abs(f(ii)) < threshold_clearance_mp);
	
}


BOOST_AUTO_TEST_CASE(patch_equality_checks)
{
	DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

	std::vector<unsigned> s1{2,3};
	std::vector<unsigned> s2{3,4};

	Patch p(s1), q(s1), r(s2);

	BOOST_CHECK_EQUAL(p,p);
	BOOST_CHECK_EQUAL(q,q);
	BOOST_CHECK_EQUAL(r,r);

	BOOST_CHECK(p!=q);
	BOOST_CHECK(p!=r);
	BOOST_CHECK(q!=r);

	BOOST_CHECK(q!=p);
	BOOST_CHECK(r!=p);
	BOOST_CHECK(r!=q);
}



BOOST_AUTO_TEST_SUITE_END() // end the patch_class test suite

