//This file is part of Bertini 2.0.
//
// python/bertini_python.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/bertini_python.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/bertini_python.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/python_common.hpp:  A common header file for all python exposure files


#ifndef Xcode_b2_python_common_hpp
#define Xcode_b2_python_common_hpp

#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python/wrapper.hpp>

#include <boost/python/operators.hpp>
#include <boost/operators.hpp>

#include <sstream>


#include <bertini2/mpfr_complex.hpp>
#include <bertini2/mpfr_extensions.hpp>
#include <bertini2/eigen_extensions.hpp>

using namespace boost::python;
typedef bertini::mpfr_float bmp;





#endif
