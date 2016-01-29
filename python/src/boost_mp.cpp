#include "boost_mp.hpp"
#include "mpfr_visitors.hpp"





namespace bertini {
	namespace python{

		using namespace boost::python;


		unsigned (*def_prec1)() = &bmp::default_precision;
		void (*def_prec2)(unsigned) = &bmp::default_precision;

		unsigned (bmp::*bmpprec1)() const= &bmp::precision;
		void (bmp::*bmpprec2)(unsigned) = &bmp::precision;

		void ExportMpfrFloat(){
			class_<bmp>("mpfr_float", init<>())
			.def(MPFRFloatVisitor<bmp>());
//			class_<bmp>("mpfr_float")
//				.def(init<bmp>())
//				.def(init<std::string>())
//				.def(init<int>())
//				.def(init<long int>())
//
//				.def(self += self)
//				.def(self + self)
//				.def(self - self)
//				.def(-self)
//				.def(self -= self)
//				.def(self * self)
//				.def(self *= self)
//				.def(self / self)
//				.def(self /= self)
//
//
//				//int
//				.def(self += other<int>())
//				.def(self + other<int>())
//				.def(other<int>() + self)
//
//				.def(self -= other<int>())
//				.def(self - other<int>())
//				.def(other<int>() - self)
//				
//				.def(other<int>() * self)
//				.def(self * other<int>())
//				.def(self *= other<int>())
//
//				.def(other<int>() / self)
//				.def(self / other<int>())
//				.def(self /= other<int>())
//
//
//
//
//
//
//
//				.def(self < self)
//				.def(self <= self)
//				.def(self > self)
//				.def(self >= self)
//
//				.def(self == self)
//
//				.def(pow(self, other<bmp>() ))
//				.def(pow(self, other<int>() ))
//
//				// .def(self_ns::str(self))
//
//				.def(abs(self))
//				
//
//				// .def("default_precision", (unsigned (bmp::*)())&bmp::default_precision) 
//				// .def("default_precision", (void (bmp::*)(unsigned int))&bmp::default_precision).staticmethod("default_precision")
//				.def("precision", bmpprec1)
//				.def("precision", bmpprec2)
//
//
//			    .def("__str__", &detail::bmp_as_str)
//			    .def("__repr__", &detail::full_precision_string)
//			    // .def("__format__", &detail::bmp_format)
//			;
//			
//			
//			
////			class_<bmp_sum>("mpfr_float_sum", no_init)
////			.def(init<bmp_sum>());
//
//			
//			
//			
//			
//			
//			def("default_precision", def_prec1);
//			def("default_precision", def_prec2);
//
//
//			def("sqrt",&detail::bmpsqrt);
//			def("abs",&detail::bmpabs);
//			def("log",&detail::bmplog);
//			def("exp",&detail::bmpexp);
//			
//			def("sin",&detail::bmpsin);
//			def("cos",&detail::bmpcos);
//			def("tan",&detail::bmptan);
//
//			def("sinh",&detail::bmpsinh);
//			def("cosh",&detail::bmpcosh);
//			def("tanh",&detail::bmptanh);
//
//			def("asin",&detail::bmpasin);
//			def("acos",&detail::bmpacos);
//			def("atan",&detail::bmpatan);

			/* These functions not yet supported on boost::multiprecision.
			def("asinh",&detail::bmpasinh);
			def("acosh",&detail::bmpacosh);
			def("atanh",&detail::bmpatanh);
			 */


			


		}

	}
}


