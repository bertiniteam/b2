//This file is part of Bertini 2.
//
//python/mpfr_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/mpfr_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/mpfr_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//
//  python/mpfr_export.cpp:  Source file for exposing all multiprecision data types, those from boost and bertini::complex.





#include "mpfr_export.hpp"



namespace bertini{
	namespace python{
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRBaseVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("__neg__",&MPFRBaseVisitor::__neg__)
			
			.def("__add__",&MPFRBaseVisitor::__add__).def("__iadd__",&MPFRBaseVisitor::__iadd__)
			.def("__sub__",&MPFRBaseVisitor::__sub__).def("__isub__",&MPFRBaseVisitor::__isub__)
			.def("__mul__",&MPFRBaseVisitor::__mul__).def("__imul__",&MPFRBaseVisitor::__imul__)
			
			.def("__str__", &MPFRBaseVisitor::__str__)
			;
			
			
			def("abs", &MPFRBaseVisitor::__abs__);
			
		}
		
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRFloatBaseVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__div__",&MPFRFloatBaseVisitor::__div__).def("__idiv__",&MPFRFloatBaseVisitor::__idiv__)
			.def("__pow__",&MPFRFloatBaseVisitor::__pow__)
			
			.def("precision", bmpprec1)
			.def("precision", bmpprec2)
			;
			
			def("exp", &MPFRFloatBaseVisitor::__exp__);
			def("log", &MPFRFloatBaseVisitor::__log__);
			def("sqrt", &MPFRFloatBaseVisitor::__sqrt__);
			
			def("sin", &MPFRFloatBaseVisitor::__sin__);
			def("cos", &MPFRFloatBaseVisitor::__cos__);
			def("tan", &MPFRFloatBaseVisitor::__tan__);
			
			def("asin", &MPFRFloatBaseVisitor::__asin__);
			def("acos", &MPFRFloatBaseVisitor::__acos__);
			def("atan", &MPFRFloatBaseVisitor::__atan__);
			
			def("sinh", &MPFRFloatBaseVisitor::__sinh__);
			def("cosh", &MPFRFloatBaseVisitor::__cosh__);
			def("tanh", &MPFRFloatBaseVisitor::__tanh__);
			
		}

		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRIntVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRIntVisitor::__add_int).def("__iadd__",&MPFRIntVisitor::__iadd_int)
			.def("__radd__",&MPFRIntVisitor::__radd_int)
			.def("__sub__",&MPFRIntVisitor::__sub_int).def("__isub__",&MPFRIntVisitor::__isub_int)
			.def("__rsub__",&MPFRIntVisitor::__rsub_int)
			.def("__mul__",&MPFRIntVisitor::__mul_int).def("__imul__",&MPFRIntVisitor::__imul_int)
			.def("__rmul__",&MPFRIntVisitor::__rmul_int)
			.def("__pow__",&MPFRIntVisitor::__pow__)
			.def("__repr__", &MPFRIntVisitor::__repr__)

			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;

		};

		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRRationalVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRRationalVisitor::__add_int).def("__iadd__",&MPFRRationalVisitor::__iadd_int)
			.def("__radd__",&MPFRRationalVisitor::__radd_int)
			.def("__sub__",&MPFRRationalVisitor::__sub_int).def("__isub__",&MPFRRationalVisitor::__isub_int)
			.def("__rsub__",&MPFRRationalVisitor::__rsub_int)
			.def("__mul__",&MPFRRationalVisitor::__mul_int).def("__imul__",&MPFRRationalVisitor::__imul_int)
			.def("__rmul__",&MPFRRationalVisitor::__rmul_int)
			.def("__div__",&MPFRRationalVisitor::__div_int).def("__idiv__",&MPFRRationalVisitor::__idiv_int)
			.def("__rdiv__",&MPFRRationalVisitor::__rdiv_int)
			
			.def("__div__",&MPFRRationalVisitor::__div__).def("__idiv__",&MPFRRationalVisitor::__idiv__)

			.def("__add__",&MPFRRationalVisitor::__add_mpint).def("__iadd__",&MPFRRationalVisitor::__iadd_mpint)
			.def("__radd__",&MPFRRationalVisitor::__radd_mpint)
			.def("__sub__",&MPFRRationalVisitor::__sub_mpint).def("__isub__",&MPFRRationalVisitor::__isub_mpint)
			.def("__rsub__",&MPFRRationalVisitor::__rsub_int)
			.def("__mul__",&MPFRRationalVisitor::__mul_mpint).def("__imul__",&MPFRRationalVisitor::__imul_mpint)
			.def("__rmul__",&MPFRRationalVisitor::__rmul_mpint)
			.def("__div__",&MPFRRationalVisitor::__div_mpint).def("__idiv__",&MPFRRationalVisitor::__idiv_mpint)

			.def("__repr__", &MPFRRationalVisitor::__repr__)
			
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;
			
		};

		
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRFloatVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRFloatBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRFloatVisitor::__add_int).def("__iadd__",&MPFRFloatVisitor::__iadd_int)
			.def("__radd__",&MPFRFloatVisitor::__radd_int)
			.def("__sub__",&MPFRFloatVisitor::__sub_int).def("__isub__",&MPFRFloatVisitor::__isub_int)
			.def("__rsub__",&MPFRFloatVisitor::__rsub_int)
			.def("__mul__",&MPFRFloatVisitor::__mul_int).def("__imul__",&MPFRFloatVisitor::__imul_int)
			.def("__rmul__",&MPFRFloatVisitor::__rmul_int)
			.def("__div__",&MPFRFloatVisitor::__div_int).def("__idiv__",&MPFRFloatVisitor::__idiv_int)
			.def("__rdiv__",&MPFRFloatVisitor::__rdiv_int)
			.def("__pow__",&MPFRFloatVisitor::__pow_int)
			
			.def("__repr__", &MPFRFloatVisitor::__repr__)
			
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)

			.def(self < int())
			.def(self <= int())
			.def(self > int())
			.def(self >= int())

			.def(int() < self)
			.def(int() <= self)
			.def(int() > self)
			.def(int() >= self)

			.def(self < double())
			.def(self <= double())
			.def(self > double())
			.def(self >= double())

			.def(double() < self)
			.def(double() <= self)
			.def(double() > self)
			.def(double() >= self)

			;
			
			
			// default_precision are defined as free functions in python
			def("default_precision", MPFRFloatVisitor<bmp>::def_prec1);
			def("default_precision", MPFRFloatVisitor<bmp>::def_prec2);

		};


		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRComplexVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRFloatBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRComplexVisitor::__add_float).def("__iadd__",&MPFRComplexVisitor::__iadd_float)
			.def("__radd__",&MPFRComplexVisitor::__radd_float)
			.def("__sub__",&MPFRComplexVisitor::__sub_float).def("__isub__",&MPFRComplexVisitor::__isub_float)
			.def("__rsub__",&MPFRComplexVisitor::__rsub_float)
			.def("__mul__",&MPFRComplexVisitor::__mul_float).def("__imul__",&MPFRComplexVisitor::__imul_float)
			.def("__rmul__",&MPFRComplexVisitor::__rmul_float)
			.def("__div__",&MPFRComplexVisitor::__div_float).def("__idiv__",&MPFRComplexVisitor::__idiv_float)
			.def("__rdiv__",&MPFRComplexVisitor::__rdiv_float)
			.def("__pow__",&MPFRComplexVisitor::__pow_int)
			.def("__pow__",&MPFRComplexVisitor::__pow_float)
			
			.def("__repr__", MPFRComplexVisitor::__repr__)
			
			.add_property("real", &MPFRComplexVisitor::get_real, &MPFRComplexVisitor::set_real)
			.add_property("imag", &MPFRComplexVisitor::get_imag, &MPFRComplexVisitor::set_imag)
			;
			
			
			// these complex-specific functions are free in python
			def("real",&real,return_value_policy<copy_const_reference>());
			def("imag",&imag,return_value_policy<copy_const_reference>());
			
			def("abs2",&MPFRBaseT::abs2);
			def("polar",&polar);
			def("norm",&MPFRBaseT::norm);
			def("conj",&MPFRBaseT::conj);
			def("arg",&arg);
			
			def("square",&square);
			def("cube",&cube);
			def("inverse", &inverse);
			def("asinh",&asinh);
			def("acosh",&acosh);
			def("atanh",&atanh);
			
		}


		void ExportMpfr()
		{
			class_<mpz_int>("mpfr_int", init<>())
			.def(init<int>())
			.def(init<mpz_int>())
			.def(MPFRIntVisitor<mpz_int>())
			;

			class_<mpq_rational>("mpfr_rational", init<>())
			.def(init<int, int>())
			.def(init<mpz_int,mpz_int>())
			.def(init<mpq_rational>())
			.def(MPFRRationalVisitor<mpq_rational>())
			;

			class_<bmp>("mpfr_float", init<>())
			.def(init<std::string>())
			.def(init<int>())
			.def(init<long int>())
			.def(init<bmp>())
			.def(MPFRFloatVisitor<bmp>())
			;
			
			class_<bertini::complex>("mpfr_complex", init<>())
			.def(init<double>())
			.def(init<mpfr_float>())
			.def(init<std::string>())
			.def(init<mpfr_float,mpfr_float>())
			.def(init<double, double>())
			.def(init<std::string, mpfr_float>())
			.def(init<mpfr_float, std::string>())
			.def(init<std::string, std::string>())
			.def(init<bertini::complex>())
			.def(MPFRComplexVisitor<bertini::complex>())
			;
			
			
			
		};

		
	} //namespace python
} // namespace bertini


