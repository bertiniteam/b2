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




// the bottom of this file is intentionally blank.  this is the 'main' .cpp file for the built boost unit test suite for bertini 2 classes.
