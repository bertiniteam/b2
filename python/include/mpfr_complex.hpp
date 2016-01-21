//This file is part of Bertini 2.0.
//
// python/mpfr_complex.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/mpfr_complex.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/mpfr_complex.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
//
//  python/mpfr_complex.hpp:  the header file for the python interface for mpfr_complex class.

#ifndef BERTINI_PYTHON_MPFR_COMPLEX_HPP
#define BERTINI_PYTHON_MPFR_COMPLEX_HPP

#include <bertini2/mpfr_complex.hpp>
#include <bertini2/function_tree/node.hpp>

#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <boost/python/wrapper.hpp>


namespace bertini{
	namespace python{
		
		void ExportMpfrComplex();

		namespace detail{

			using cplx = bertini::complex;

			
			inline
			cplx cplxsqrt(cplx const& n)
			{
				return sqrt(n);
			}


			// a free function getting a number as a string
			inline
			std::string cplx_as_str(cplx const& n)
			{
				std::stringstream s;
				s << n;
				return s.str();
			}

			// a free function getting a number as a string, to full precision
			inline
			std::string cplx_full_precision_string(cplx const& n)
			{	
				std::stringstream s;
				s << "(" << real(n).str(0,std::ios::scientific) << ", " << imag(n).str(0,std::ios::scientific) << ")";
				return s.str();
			}


		} // namespace detail
	} // namespace python
} // namespace bertini
#endif


