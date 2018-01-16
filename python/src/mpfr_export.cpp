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
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Danielle Brake
//  University of Wisconsin - Eau Claire
//  Fall 2017, Spring 2018
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
		
		template<typename T>
		template<typename PyClass>
		void PrecisionVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def("precision", get_prec)
			.def("precision", set_prec)
			;
		}

		template<typename T>
		template<typename PyClass>
		void RealStrVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def("__str__", &RealStrVisitor::__str__)
			.def("__repr__", &RealStrVisitor::__repr__)
			;
		}

		template<typename T>
		template<typename PyClass>
		void EqualitySelfVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def(self == self)
			.def(self != self)
			;
		}

		template<typename T, typename S>
		template<typename PyClass>
		void EqualityVisitor<T,S>::visit(PyClass& cl) const
		{
			cl
			.def(self == other<S>())
			.def(self != other<S>())

			.def(other<S>() == self)
			.def(other<S>() != self)
			;
		}


		template<typename T, typename S>
		template<typename PyClass>
		void RingVisitor<T,S>::visit(PyClass& cl) const
		{
			cl
			.def("__add__",&RingVisitor::__add__)
			.def("__iadd__",&RingVisitor::__iadd__)
			.def("__radd__",&RingVisitor::__radd__)

			.def("__sub__",&RingVisitor::__sub__)
			.def("__isub__",&RingVisitor::__isub__)
			.def("__rsub__",&RingVisitor::__rsub__)

			.def("__mul__",&RingVisitor::__mul__)
			.def("__imul__",&RingVisitor::__imul__)
			.def("__rmul__",&RingVisitor::__rmul__)
			;
		}


		template<typename T>
		template<typename PyClass>
		void RingSelfVisitor<T>::visit(PyClass& cl) const
		{
			cl

			.def("__add__",&RingSelfVisitor::__add__)
			.def("__iadd__",&RingSelfVisitor::__iadd__)

			.def("__sub__",&RingSelfVisitor::__sub__)
			.def("__isub__",&RingSelfVisitor::__isub__)

			.def("__mul__",&RingSelfVisitor::__mul__)
			.def("__imul__",&RingSelfVisitor::__imul__)

			.def("__neg__",&RingSelfVisitor::__neg__)
			;

			
		}


		template<typename T>
		template<typename PyClass>
		void RealFreeVisitor<T>::visit(PyClass& cl) const
		{
			def("abs", &RealFreeVisitor::__abs__); // free
		}



		template<typename T, typename S>
		template<typename PyClass>
		void FieldVisitor<T,S>::visit(PyClass& cl) const
		{
			cl
			.def("__div__",&FieldVisitor::__div__)
			.def("__idiv__",&FieldVisitor::__idiv__)
			.def("__rdiv__",&FieldVisitor::__rdiv__)

			.def("__truediv__",&FieldVisitor::__div__)
			.def("__itruediv__",&FieldVisitor::__idiv__)
			.def("__rtruediv__",&FieldVisitor::__rdiv__)
			.def(RingVisitor<T,S>())
			;
		}


		template<typename T>
		template<typename PyClass>
		void FieldSelfVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def("__div__",&FieldSelfVisitor::div)
			.def("__idiv__",&FieldSelfVisitor::idiv)

			.def("__truediv__",&FieldSelfVisitor::div)
			.def("__itruediv__",&FieldSelfVisitor::idiv)

			.def(RingSelfVisitor<T>())
			;
		}

		template<typename T, typename S>
		template<typename PyClass>
		void PowVisitor<T,S>::visit(PyClass& cl) const
		{
			cl
			.def("__pow__",&PowVisitor::__pow__)
			;
		}

		template<typename T, typename S>
		template<typename PyClass>
		void GreatLessVisitor<T,S>::visit(PyClass& cl) const
		{
			cl
			.def(self < other<S>())
			.def(self <= other<S>())
			.def(self > other<S>())
			.def(self >= other<S>())

			.def(other<S>() < self)
			.def(other<S>() <= self)
			.def(other<S>() > self)
			.def(other<S>() >= self)
			;
		}

		template<typename T>
		template<typename PyClass>
		void GreatLessSelfVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			;
		}


		template<typename T>
		template<typename PyClass>
		void TranscendentalVisitor<T>::visit(PyClass& cl) const
		{
			def("exp", &TranscendentalVisitor::__exp__);
			def("log", &TranscendentalVisitor::__log__);
			def("sqrt", &TranscendentalVisitor::__sqrt__);
			
			def("sin", &TranscendentalVisitor::__sin__);
			def("cos", &TranscendentalVisitor::__cos__);
			def("tan", &TranscendentalVisitor::__tan__);
			
			def("asin", &TranscendentalVisitor::__asin__);
			def("acos", &TranscendentalVisitor::__acos__);
			def("atan", &TranscendentalVisitor::__atan__);
			
			def("sinh", &TranscendentalVisitor::__sinh__);
			def("cosh", &TranscendentalVisitor::__cosh__);
			def("tanh", &TranscendentalVisitor::__tanh__);

			def("asinh",&TranscendentalVisitor::__asinh__);
			def("acosh",&TranscendentalVisitor::__acosh__);
			def("atanh",&TranscendentalVisitor::__atanh__);
		}


		template<typename T>
		template<class PyClass>
		void ComplexVisitor<T>::visit(PyClass& cl) const
		{
			// MPFRFloatBaseVisitor<T>().visit(cl);
			
			cl			
			.add_property("real", &ComplexVisitor::get_real, &ComplexVisitor::set_real)
			.add_property("imag", &ComplexVisitor::get_imag, &ComplexVisitor::set_imag)
			
			.def("__str__", &ComplexVisitor::__str__)
			.def("__repr__", &ComplexVisitor::__repr__)
			;
			
			
			// these complex-specific functions are free in python
			def("real",&real,return_value_policy<copy_const_reference>());
			def("imag",&imag,return_value_policy<copy_const_reference>());
			
			// and then a few more free functions
			def("abs2",&T::abs2);
			def("polar",&polar);
			def("norm",&T::norm);
			def("conj",&T::conj);
			def("arg",&arg);
			
			def("square",&square);
			def("cube",&cube);
			def("inverse", &inverse);

			def("abs", &ComplexVisitor::__abs__); // free
		}










		void ExposeFreeNumFns()
		{
			unsigned (*def_prec1)() = &bertini::DefaultPrecision;
			void (*def_prec2)(unsigned) = &bertini::DefaultPrecision;

			def("default_precision", def_prec1);
			def("default_precision", def_prec2);
		}





		void ExposeInt()
		{
			using T = mpz_int;

			class_<mpz_int>("Int", init<>())
			.def(init<int>())
			.def(init<T>())
			.def(RealStrVisitor<T>())
			.def(RingSelfVisitor<T>())
			.def(PowVisitor<T,int>())
			.def(GreatLessSelfVisitor<T>())
			.def(GreatLessVisitor<T,int>())

			.def(EqualitySelfVisitor<T>())
			.def(EqualityVisitor<T, int>())

			.def(RealFreeVisitor<T>())
			;
		}


		

		
		void ExposeRational()
		{
			using T = mpq_rational;

			class_<mpq_rational>("Rational", init<>())
			.def(init<int>())
			.def(init<int, int>())
			.def(init<mpz_int>())
			.def(init<mpz_int,mpz_int>())
			.def(init<mpq_rational>())
			.def(RealStrVisitor<T>())
			.def(FieldSelfVisitor<T>())
			.def(FieldVisitor<T, mpz_int>())
			// .def(PowVisitor<T,int>()) // deliberately commented out... 
										 // pow(Q,Z) not defined...
			.def(GreatLessSelfVisitor<T>())
			.def(GreatLessVisitor<T,int>())
			.def(GreatLessVisitor<T,mpz_int>())

			.def(EqualitySelfVisitor<T>())
			.def(EqualityVisitor<T, int>())
			.def(EqualityVisitor<T, mpz_int>())

			.def(RealFreeVisitor<T>())
			;
		}




		void ExposeFloat()
		{
			using T = bmp;

			class_<T>("Float", init<>())
			.def(init<std::string>())
			.def(init<int>())
			.def(init<long int>())
			.def(init<T>())

			.def(init<mpz_int>())

			.def(RealStrVisitor<T>())
			.def(PrecisionVisitor<T>())

			.def(FieldSelfVisitor<T>())
			
			.def(FieldVisitor<T, int>())
			.def(FieldVisitor<T, mpz_int>())
			.def(FieldVisitor<T, mpq_rational>())
			
			.def(PowVisitor<T,T>())
			.def(PowVisitor<T,int>())
			.def(TranscendentalVisitor<T>())

			.def(GreatLessSelfVisitor<T>())
			.def(GreatLessVisitor<T,int>())
			.def(GreatLessVisitor<T,double>())

			.def(RealFreeVisitor<T>())
			;
		}




		void ExposeComplex()
		{

			using T = bertini::complex;

			class_<T>("Complex", init<>())
			.def(init<double>())
			.def(init<mpfr_float>())
			.def(init<std::string>())
			.def(init<mpfr_float,mpfr_float>())
			.def(init<double, double>())
			.def(init<std::string, mpfr_float>())
			.def(init<mpfr_float, std::string>())
			.def(init<std::string, std::string>())
			.def(init<T>())

			.def(init<mpz_int>())
			.def(init<mpz_int, mpz_int>())
			
			.def(ComplexVisitor<T>())

			.def(FieldSelfVisitor<T>())

			.def(FieldVisitor<T, mpz_int>())
			.def(FieldVisitor<T, mpq_rational>())
			.def(FieldVisitor<T, mpfr_float>())

			.def(FieldVisitor<T, int>())

			.def(PowVisitor<T,T>())
			.def(PowVisitor<T,int>())
			.def(PowVisitor<T,mpfr_float>())

			.def(TranscendentalVisitor<T>())

			.def(PrecisionVisitor<T>())
			;
		}

		void ExportMpfr()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".multiprec");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("multiprec") = new_submodule;
			scope new_submodule_scope = new_submodule;

			ExposeInt();
			ExposeFloat();
			ExposeRational();
			ExposeComplex();	

			ExposeFreeNumFns();	
		};

		
	} //namespace python
} // namespace bertini


