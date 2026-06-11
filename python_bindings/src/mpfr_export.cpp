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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin - Eau Claire
//  Fall 2017, Spring 2018, Spring 2026
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
			cl.add_property("precision",
				get_prec, set_prec,
				"get/set the precision of this variable-precision number, in digits.  remember, the system knows not where your number came from, so upsampling will NOT add more correct digits.");
		}

		template<typename T>
		template<typename PyClass>
		void RealStrVisitor<T>::visit(PyClass& cl) const
		{
			cl
			.def("__str__", &RealStrVisitor::__str__, (arg("self")))
			.def("__repr__", &RealStrVisitor::__repr__, (arg("self")))
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
			.def("__add__",&RingVisitor::__add__, "addition")
			.def("__iadd__",&RingVisitor::__iadd__, "addition")
			.def("__radd__",&RingVisitor::__radd__, "addition")

			.def("__sub__",&RingVisitor::__sub__, "subtraction")
			.def("__isub__",&RingVisitor::__isub__, "subtraction")
			.def("__rsub__",&RingVisitor::__rsub__, "subtraction")

			.def("__mul__",&RingVisitor::__mul__, "multiplication")
			.def("__imul__",&RingVisitor::__imul__, "multiplication")
			.def("__rmul__",&RingVisitor::__rmul__, "multiplication")
			;
		}


		template<typename T>
		template<typename PyClass>
		void RingSelfVisitor<T>::visit(PyClass& cl) const
		{
			cl

			.def("__add__",&RingSelfVisitor::__add__, "addition")
			.def("__iadd__",&RingSelfVisitor::__iadd__, "in-place addition")

			.def("__sub__",&RingSelfVisitor::__sub__, "subtraction")
			.def("__isub__",&RingSelfVisitor::__isub__, "in-place subtraction")

			.def("__mul__",&RingSelfVisitor::__mul__, "multiplication")
			.def("__imul__",&RingSelfVisitor::__imul__, "in-place multiplication")

			.def("__neg__",&RingSelfVisitor::__neg__, "negation")
			;


		}


		template<typename T>
		template<typename PyClass>
		void RealFreeVisitor<T>::visit(PyClass& cl) const
		{
			def("abs", &RealFreeVisitor::__abs__, (arg("val")), "absolute value"); // free
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
			.def("__div__",&FieldSelfVisitor::div, "division")
			.def("__idiv__",&FieldSelfVisitor::idiv, "division")

			.def("__truediv__",&FieldSelfVisitor::div, "division")
			.def("__itruediv__",&FieldSelfVisitor::idiv, "division")

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
			def("exp", &TranscendentalVisitor::__exp__, (arg("val")), "exponential, base e");
			def("log", &TranscendentalVisitor::__log__, (arg("val")), "natural log");
			def("sqrt", &TranscendentalVisitor::__sqrt__, (arg("val")), "square root");

			def("sin", &TranscendentalVisitor::__sin__, (arg("val")), "sine");
			def("cos", &TranscendentalVisitor::__cos__, (arg("val")), "cosine");
			def("tan", &TranscendentalVisitor::__tan__, (arg("val")), "tangent");

			def("asin", &TranscendentalVisitor::__asin__, (arg("val")), "arcsine");
			def("acos", &TranscendentalVisitor::__acos__, (arg("val")), "arccosine");
			def("atan", &TranscendentalVisitor::__atan__, (arg("val")), "arctangent");

			def("sinh", &TranscendentalVisitor::__sinh__, (arg("val")), "hyperbolic sine");
			def("cosh", &TranscendentalVisitor::__cosh__, (arg("val")), "hyperbolic cosine");
			def("tanh", &TranscendentalVisitor::__tanh__, (arg("val")), "hyperbolic tangent");

			def("asinh",&TranscendentalVisitor::__asinh__, (arg("val")), "hyperbolic arcsine");
			def("acosh",&TranscendentalVisitor::__acosh__, (arg("val")), "hyperbolic arccosine");
			def("atanh",&TranscendentalVisitor::__atanh__, (arg("val")), "hyperbolic arctangent");
		}


		template<typename T>
		template<class PyClass>
		void ComplexVisitor<T>::visit(PyClass& cl) const
		{
			// MPFRFloatBaseVisitor<T>().visit(cl);

			cl
			.add_property("real", &ComplexVisitor::get_real, &ComplexVisitor::set_real,"the real part of the complex number")
			.add_property("imag", &ComplexVisitor::get_imag, &ComplexVisitor::set_imag,"the imaginary part of the complex number")

			.def("__str__", &ComplexVisitor::__str__, (arg("self")), "convert to string")
			.def("__repr__", &ComplexVisitor::__repr__, (arg("self")), "convert to string")
			;


			// these complex-specific functions are free in python
			using boost::multiprecision::real;
			using boost::multiprecision::imag;

			mpfr_float (*reeeal)(const T&) = &boost::multiprecision::real;
			mpfr_float (*imaaag)(const T&) = &boost::multiprecision::real;
			def("real",reeeal, (arg("val")), "get the real part"); //,return_value_policy<copy_const_reference>()
			def("imag",imaaag, (arg("val")), "get the imaginary part"); //,return_value_policy<copy_const_reference>()

			// and then a few more free functions
			// def("abs2",&T::abs2);

			mpfr_complex (*pooolar)(const mpfr_float&,const mpfr_float&) = &boost::multiprecision::polar;

			def("polar",pooolar, "construct from polar form");
			// def("norm",&T::norm);

			T (*conjjj)(const T&) = +[](const T& x) -> T { return boost::multiprecision::conj(x); };
			def("conj",conjjj, "complex conjugate");

			mpfr_float (*aaaarg)(const T&) = &boost::multiprecision::arg;
			def("arg",aaaarg, "the argument, or the angle from 0.  beware the branch cut.");

			// def("square",&square);
			// def("cube",&cube);
			// def("inverse", &inverse);

			def("abs", &ComplexVisitor::__abs__, "the magnitude of a complex number"); // free
		}



		template <typename T>
		unsigned get_precision_vector(Eigen::Ref<Vec<T>> x)
		{
			return bertini::Precision(x);
		}


#define IMPLICITLY_CONVERTIBLE(T1,T2) \
  boost::python::implicitly_convertible<T1,T2>();



		void ExposeFreeNumFns()
		{
			unsigned (*def_prec1)() = &bertini::DefaultPrecision;
			void (*def_prec2)(unsigned) = &bertini::DefaultPrecision;

			// Seed thread_default_precision on the interpreter thread so that
			// boost::python const& extractors don't construct mpfr temporaries
			// at precision 0 on Boost >= 1.87 (which aborts inside mpfr_init2).
			// Read the static default directly; fall back to 50 (Boost's library
			// default) in case the getter itself delegates to the unset thread-local.
			{
				auto p = mpfr_float::default_precision();
				if (p == 0) p = 50;
				mpfr_float::thread_default_precision(p);
				mpfr_complex::thread_default_precision(p);
			}

			def("default_precision", def_prec1, "get the default precision for variable-precision numbers.  is digits, not bits.");
			def("default_precision", def_prec2, "set the default precision for variable-precision numbers.  should be a positive number.  is digits, not bits.");
		}





		void ExposeInt()
		{
			using T = mpz_int;

			class_<mpz_int>("Int", init<>("Default Construct an arbitrary-precision integer"))
			.def(init<int>((arg("self"),arg("val")),"Construct an arbitrary-precision integer from an integer."))
			.def(init<T>((arg("self"),arg("val")),"Construct an arbitrary-precision integer from another."))
			.def(init<std::string>((arg("self"),arg("val")),"Construct an arbitrary-precision integer from a string of digits."))
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

			class_<mpq_rational>("Rational", init<>("Default Construct an arbitrary-precision rational number"))
			.def(init<int>((arg("self"),arg("val")),"Construct an arbitrary-precision rational number from an integer."))
			.def(init<int, int>((arg("self"),arg("numerator"), arg("denominator")),"Construct an arbitrary-precision rational number from a pair of integers."))
			.def(init<mpz_int>((arg("self"),arg("val")),"Construct an arbitrary-precision rational number from an arbitrary-precision integer."))
			.def(init<mpz_int,mpz_int>((arg("self"),arg("numerator"),arg("denominator")),"Construct an arbitrary-precision rational number from a pair of arbitrary-precision integers."))
			.def(init<std::string>((arg("self"),arg("val")),"Construct an arbitrary-precision rational number from a string, e.g. '1/3'."))
			.def(init<mpq_rational>((arg("self"),arg("val")),"Construct an arbitrary-precision rational number from an arbitrary-precision integer."))
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
			using T = mpfr_float;

			class_<T>("Float", init<>("Default Construct a variable-precision float"))
			.def(init<std::string>((arg("self"),arg("val")),"Construct a variable-precision float from a string.  The best way."))
			.def(init<long int>((arg("self"),arg("val")),"Construct a variable-precision float from a regular old integer."))
			.def(init<T>((arg("self"),arg("val")),"Construct a variable-precision float from another."))

			.def(init<mpz_int>((arg("self"),arg("val")),"Construct an variable-precision float from an arbitrary-precision integer."))

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

			.def(EqualitySelfVisitor<T>())
			.def(EqualityVisitor<T, int>())
			.def(EqualityVisitor<T, mpz_int>())

			.def(RealFreeVisitor<T>())
			;


			eigenpy::registerNewType<T>();
			eigenpy::HardenSetitem<T>(); // zero slots before assignment — see eigenpy_interaction.hpp & ADR-0003
			eigenpy::HardenDotfunc<T>(); // np.dot/np.inner guard — see eigenpy_interaction.hpp
			// guarded loops (real type — orderings included); eigenpy's registerCommonUfunc
			// loops read input slots unguarded and crash on never-written np.zeros/np.empty slots.
			eigenpy::registerGuardedUfunct<T, true>();

			// you can convert from integer types with no fear
			eigenpy::registerCast<long,T>(true);
			eigenpy::registerCast<int,T>(true);
			eigenpy::registerCast<int64_t,T>(true);

			// but you can never convert TO integer types.  so these are commented out.
			// eigenpy::registerCast<T,long>(true);
			// eigenpy::registerCast<T,int>(true);
			// eigenpy::registerCast<T,int64_t>(true);

			// unsafe.  so the argument is false.
			// you can ask for the conversions, but you probably shouldn't.
			// both directions are scary.
			eigenpy::registerCast<T,double>(false);
			eigenpy::registerCast<double,T>(false);


			IMPLICITLY_CONVERTIBLE(int,T);
			IMPLICITLY_CONVERTIBLE(long,T);
			IMPLICITLY_CONVERTIBLE(int64_t,T);

			// do not allow implicit conversion, because it is a potential source of problems.
			// because 0.1 as a float64 does NOT convert to 0.1 as a variable precision number.
			// the user should use strings to guarantee matching.
			// that is, leave commented-out.  silviana, 2026.04.14
			// IMPLICITLY_CONVERTIBLE(double,T);


			eigenpy::EigenToPyConverter<Vec<T>>::registration();
			eigenpy::EigenToPyConverter<Mat<T>>::registration();
			eigenpy::EigenFromPyConverter<Vec<T>>::registration();
			eigenpy::EigenFromPyConverter<Mat<T>>::registration();
		}

		size_t get_default_align(){return EIGENPY_DEFAULT_ALIGN_BYTES;}


		void ExposeComplex()
		{

			using T = bertini::mpfr_complex;

			class_<T>("Complex", init<>())
			.def(init<double>((arg("self"),arg("real")),"Construct variable-precision complex number from a double, with 0 imaginary part. do this with caution, as 0.1 is not what you think it is -- there's noise at the end.")) // this should probably be made an explicit constructor rather than implicit
			.def(init<mpfr_float>((arg("self"),arg("real")),"Construct variable-precision complex number from a variable-precision float, with 0 imaginary part"))
			.def(init<std::string>((arg("self"),arg("real")),"Construct variable-precision complex number from a string, with 0 imaginary part"))
			.def(init<mpfr_float,mpfr_float>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a pair of variable-precision floats"))
			.def(init<double, double>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a pair of doubles.  do this with caution, as 0.1 is not what you think it is -- there's noise at the end.")) // this should probably be made an explicit constructor rather than implicit
			.def(init<std::string, mpfr_float>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a string and a variable-precision float"))
			.def(init<mpfr_float, std::string>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a variable-precision float and a string"))
			.def(init<std::string, std::string>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a pair of strings.  the best way to construct one and be sure you have padded with zeros to the end, in the current working precision"))

			.def(init<T>((arg("self"),arg("value")),"Construct variable-precision complex number from another one"))

			.def(init<mpz_int>((arg("self"),arg("real")),"Construct variable-precision complex number from an arbitrary-precision integer, with 0 imaginary part"))
			.def(init<mpz_int, mpz_int>((arg("self"),arg("real"),arg("imag")),"Construct variable-precision complex number from a pair of arbitrary-precision integers"))

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

			.def(EqualitySelfVisitor<T>())
			.def(EqualityVisitor<T, mpfr_float>())
			.def(EqualityVisitor<T, int>())
			;


			eigenpy::registerNewType<T>();
			eigenpy::HardenSetitem<T>(); // zero slots before assignment — see eigenpy_interaction.hpp & ADR-0003
			eigenpy::HardenDotfunc<T>(); // np.dot/np.inner guard — see eigenpy_interaction.hpp
			eigenpy::registerUfunct_without_comparitors<T>();


			// you can safely convert from integer types to Complex's, there's no loss possible
			eigenpy::registerCast<long,T>(true);
			eigenpy::registerCast<int,T>(true);
			eigenpy::registerCast<int64_t,T>(true);

			// you can never convert TO integer types from Complex, so these are commented out.
			// do not comment them in.
			// eigenpy::registerCast<T,long>(false);
			// eigenpy::registerCast<T,int>(false);
			// eigenpy::registerCast<T,int64_t>(false);

			// these conversions are unsafe.  you can, but you probably shouldn't
			eigenpy::registerCast<T,double>(false);
			eigenpy::registerCast<double,T>(false);

			// it's ok to convert from variable precision Float to Complex, that's ok!
			eigenpy::registerCast<mpfr_float,T>(true);

			IMPLICITLY_CONVERTIBLE(int,T);
			IMPLICITLY_CONVERTIBLE(long,T);
			IMPLICITLY_CONVERTIBLE(int64_t,T);

			// this is a general python conversion, not an eigenpy conversion.  it's ok to convert from reals to complexes.
			IMPLICITLY_CONVERTIBLE(mpfr_float,T);

			// BUT!!
			// do not allow implicit conversion, because it is a potential source of problems.
			// because 0.1 as a float64 does NOT convert to 0.1 as a variable precision number.
			// the user should use strings to guarantee matching.
			// that is, leave commented-out.  silviana, 2026.04.14
			// IMPLICITLY_CONVERTIBLE(double,T);

			eigenpy::EigenToPyConverter<Vec<T>>::registration();
			eigenpy::EigenToPyConverter<Mat<T>>::registration();
			eigenpy::EigenFromPyConverter<Vec<T>>::registration();
			eigenpy::EigenFromPyConverter<Mat<T>>::registration();

			eigenpy::exposeType<T>();
			eigenpy::exposeType<T, Eigen::RowMajor>();

			boost::python::def("precision", &get_precision_vector<mpfr_complex>, "get the precision of a vector of complexes");

			boost::python::def("default_align_bytes", &get_default_align);

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


#undef IMPLICITLY_CONVERTIBLE

	} //namespace python
} // namespace bertini


