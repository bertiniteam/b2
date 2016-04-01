//This file is part of Bertini 2.0.
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
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
// class_test.cpp:  main source file for class testing executable for Bertini2






// the purpose of this file is the three non-comment lines.  it has no other purpose than to provide a place for them.



//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Class Testing"
#include <boost/test/unit_test.hpp>


#include "bertini2/num_traits.hpp"

using dbl = bertini::dbl;
using mpfr = bertini::complex;

double relaxed_threshold_clearance_d = 1e-14;
double threshold_clearance_d = 1e-15;

unsigned CLASS_TEST_MPFR_DEFAULT_DIGITS = 30;
bertini::mpfr_float threshold_clearance_mp("1e-27");

std::string xstr_real = "3.1";
std::string xstr_imag = "4.1";
std::string ystr_real = "8.8";
std::string ystr_imag = "9.9";
std::string zstr_real = "0.8";
std::string zstr_imag = "-1.7";

std::string astr_real = "3.4";
std::string astr_imag = "5.6";
std::string bstr_real = "-0.2";
std::string bstr_imag = "-2.1";
std::string pstr_real = "0.8";
std::string pstr_imag = "-1.7";



dbl xnum_dbl(std::stod(xstr_real), std::stod(xstr_imag));
dbl ynum_dbl(std::stod(ystr_real), std::stod(ystr_imag));
dbl znum_dbl(std::stod(zstr_real), std::stod(zstr_imag));
dbl anum_dbl(std::stod(astr_real), std::stod(astr_imag));
dbl bnum_dbl(std::stod(bstr_real), std::stod(bstr_imag));
dbl pnum_dbl(std::stod(pstr_real), std::stod(pstr_imag));

mpfr xnum_mpfr(xstr_real, xstr_imag);
mpfr ynum_mpfr(ystr_real, ystr_imag);
mpfr znum_mpfr(zstr_real, zstr_imag);
mpfr anum_mpfr(astr_real, astr_imag);
mpfr bnum_mpfr(bstr_real, bstr_imag);
mpfr pnum_mpfr(pstr_real, pstr_imag);

// the bottom of this file is intentionally blank.  this is the 'main' .cpp file for the built boost unit test suite for bertini 2 classes.
