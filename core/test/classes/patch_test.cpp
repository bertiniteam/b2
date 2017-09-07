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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

//  patch_test.cpp
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
// also modified by
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

/**
\file patch_test.cpp Unit testing for the bertini::Patch class.
*/


#include <boost/test/unit_test.hpp>
#include "bertini2/system/patch.hpp"




#include "externs.hpp"




BOOST_AUTO_TEST_SUITE(patch_class)


using bertini::DefaultPrecision;

template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;
using Patch = bertini::Patch;

using dbl = bertini::dbl;
using mpfr = bertini::mpfr;

using mpfr_float = bertini::mpfr_float;


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

	Vec<dbl> v(5);
	v << dbl(1),  dbl(1),  dbl(1),  dbl(1),  dbl(1);

	p.Precision(16);
	
	auto f = p.Eval(v);

	BOOST_CHECK_EQUAL(f.size(),2);
}


BOOST_AUTO_TEST_CASE(patch_jacobian_two_variable_groups_prec16)
{
	std::vector<unsigned> s{2,3};

	Patch p(s);

	Vec<dbl> v(5);
	v << dbl(1),  dbl(1),  dbl(1),  dbl(1),  dbl(1);

	p.Precision(16);

	auto J = p.Jacobian(v);
	BOOST_CHECK_EQUAL(J.rows(),2);
	BOOST_CHECK_EQUAL(J.cols(),5);

	BOOST_CHECK_EQUAL(J(0,2),dbl(0));
	BOOST_CHECK_EQUAL(J(0,3),dbl(0));
	BOOST_CHECK_EQUAL(J(0,4),dbl(0));

	BOOST_CHECK_EQUAL(J(1,0),dbl(0));
	BOOST_CHECK_EQUAL(J(1,1),dbl(0));
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

	Vec<dbl> v(5);
	v << dbl(1),  dbl(1),  dbl(1),  dbl(1),  dbl(1);

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

