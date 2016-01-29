#include "mpfr_complex.hpp"


using namespace boost::python;

using mpfr_float = bertini::mpfr_float;



BOOST_PYTHON_FUNCTION_OVERLOADS(MpfrComplexSetRealOverload, set_real, 1, 1);
BOOST_PYTHON_FUNCTION_OVERLOADS(MpfrComplexSetImaOverload, set_imag, 1, 1);

namespace bertini{
	namespace python{


		

		

		// typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0> , boost::multiprecision::et_off> bmp;

		unsigned (bertini::complex::*cplxprec1)() const= &bertini::complex::precision;
		void (bertini::complex::*cplxprec2)(unsigned) = &bertini::complex::precision;


		void ExportMpfrComplex()
		{
			using complex = bertini::complex;

			class_<complex>("mpfr_complex")
				.def(init<double>())
				.def(init<mpfr_float>())
				.def(init<std::string>())
				.def(init<mpfr_float,mpfr_float>())
				.def(init<double, double>())
				.def(init<std::string, mpfr_float>())
				.def(init<mpfr_float, std::string>())
				.def(init<std::string, std::string>())


				// homogeneous arithmetic operators

				// addition
				.def(self += self)
				.def(self + self)

				// subtraction
				.def(self -= self)
				.def(self - self)
				
				// multiplication
				.def(self *= self)
				.def(self * self)

				// division
				.def(self /= self)
				.def(self / self)

				// negation
				.def(-self)

				///////
				// and now the mixed type arithmetic definitions
				////////

				//mpfr_float
				.def(self += other<mpfr_float>())
				.def(self + other<mpfr_float>())
				.def(other<mpfr_float>() + self)

				.def(self -= other<mpfr_float>())
				.def(self - other<mpfr_float>())
				.def(other<mpfr_float>() - self)
				
				.def(other<mpfr_float>() * self)
				.def(self * other<mpfr_float>())
				.def(self *= other<mpfr_float>())

				.def(other<mpfr_float>() / self)
				.def(self / other<mpfr_float>())
				.def(self /= other<mpfr_float>())

			
				


				// powers
				.def(pow(self, self ))
				.def(pow(self, other<mpfr_float>() ))
				.def(pow(self, other<int>() ))

				// .def(self_ns::str(self))

				.def(abs(self))
				

				// comparitors
				.def(self == self)

				// .def("default_precision", (unsigned (bmp::*)())&bmp::default_precision) 
				// .def("default_precision", (void (bmp::*)(unsigned int))&bmp::default_precision).staticmethod("default_precision")
				.def("precision", cplxprec1)
				.def("precision", cplxprec2)


			    .def("__str__", &detail::cplx_as_str)
			    .def("__repr__", &detail::cplx_full_precision_string)


	
				.add_property("real", &get_real, &set_real)
				.add_property("imag", &get_imag, &set_imag)



				//MpfrComplexSetRealOverload(args("new_real", "docstring"))
				//, MpfrComplexRealOverload(args("new_real"), "real's docstring")
				// .def("reset", &bertini::node::Node::Reset, &bertini::python::NodeWrap::default_Reset)

				// .def("Eval", &bertini::node::Node::Eval<dbl>,
				// 	NodeEvalOverloadsDbl(args("DiffVar"), "Eval's docstring"))
			;


			def("sqrt",&detail::cplxsqrt);

			def("real",&real);
			def("imag",&imag);

			def("abs",&abs);
			def("abs2",&abs2);
			def("sqrt",&sqrt);

			def("polar",&polar);
			def("norm",&complex::norm);
			def("conj",&conj);
			def("arg",&arg);
			def("inverse",&inverse);
			def("square",&square);
			def("cube",&cube);
			def("polar",&polar);
			def("exp",&exp);
			def("sin",&sin);
			def("cos",&cos);
			def("tan",&tan);

			def("sinh",&sinh);
			def("cosh",&cosh);
			def("tanh",&tanh);

			def("asin",&asin);
			def("acos",&acos);
			def("atan",&atan);

			def("asinh",&asinh);
			def("acosh",&acosh);
			def("atanh",&atanh);


			def("log",&log);


		}

	}
}
