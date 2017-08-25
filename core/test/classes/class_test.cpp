//This file is part of Bertini 2.
//
//class_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//class_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with class_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

// class_test.cpp:  main source file for class testing executable for Bertini2






// the purpose of this file is the three non-comment lines.  it has no other purpose than to provide a place for them.



//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Class Testing"
#include <boost/test/unit_test.hpp>


#include "bertini2/num_traits.hpp"

#include "externs.hpp"

using dbl = bertini::dbl;
using mpfr = bertini::complex;

const double relaxed_threshold_clearance_d = 1e-14;
const double threshold_clearance_d = 1e-15;

unsigned const CLASS_TEST_MPFR_DEFAULT_DIGITS = 50;
bertini::mpfr_float threshold_clearance_mp("1e-27");

const std::string xstr_real = "3.1";
const std::string xstr_imag = "4.1";
const std::string ystr_real = "8.8";
const std::string ystr_imag = "9.9";
const std::string zstr_real = "0.8";
const std::string zstr_imag = "-1.7";

const std::string astr_real = "3.4";
const std::string astr_imag = "5.6";
const std::string bstr_real = "-0.2";
const std::string bstr_imag = "-2.1";
const std::string pstr_real = "0.8";
const std::string pstr_imag = "-1.7";



const dbl xnum_dbl(std::stod(xstr_real), std::stod(xstr_imag));
const dbl ynum_dbl(std::stod(ystr_real), std::stod(ystr_imag));
const dbl znum_dbl(std::stod(zstr_real), std::stod(zstr_imag));
const dbl anum_dbl(std::stod(astr_real), std::stod(astr_imag));
const dbl bnum_dbl(std::stod(bstr_real), std::stod(bstr_imag));
const dbl pnum_dbl(std::stod(pstr_real), std::stod(pstr_imag));

mpfr xnum_mpfr;
mpfr ynum_mpfr;
mpfr znum_mpfr;
mpfr anum_mpfr;
mpfr bnum_mpfr;
mpfr pnum_mpfr;

/**
This quick struct acts only to initialize the mpfr's to the correct precision, with correct values.  Brent discovered on his machine that these variables were in the incorrect precision, and it caused many tests to fail.  Hence, this fixture.
*/
struct ClassTestInit{
	ClassTestInit()
	{
		bertini::DefaultPrecision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		xnum_mpfr.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		ynum_mpfr.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		znum_mpfr.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		anum_mpfr.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		bnum_mpfr.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);
		threshold_clearance_mp.precision(CLASS_TEST_MPFR_DEFAULT_DIGITS);

		xnum_mpfr = mpfr(xstr_real, xstr_imag);
		ynum_mpfr = mpfr(ystr_real, ystr_imag);
		znum_mpfr = mpfr(zstr_real, zstr_imag);
		anum_mpfr = mpfr(astr_real, astr_imag);
		bnum_mpfr = mpfr(bstr_real, bstr_imag);
		pnum_mpfr = mpfr(pstr_real, pstr_imag);

		threshold_clearance_mp = bertini::mpfr_float("1e-27");
	}
};


BOOST_GLOBAL_FIXTURE( ClassTestInit );


// the bottom of this file is intentionally blank.  this is the 'main' .cpp file for the built boost unit test suite for bertini 2 classes.
