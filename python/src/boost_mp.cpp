#include "boost_mp.hpp"





namespace bertini {
	namespace python{

		using namespace boost::python;


		unsigned (*def_prec1)() = &bmp::default_precision;
		void (*def_prec2)(unsigned) = &bmp::default_precision;

		unsigned (bmp::*bmpprec1)() const= &bmp::precision;
		void (bmp::*bmpprec2)(unsigned) = &bmp::precision;

		void ExportMpfrFloat(){
			class_<bmp>("mpfr_float")
				.def(init<std::string>())
				.def(init<int>())
				.def(init<long int>())

				.def(self += self)
				.def(self + self)
				.def(self - self)
				.def(-self)
				.def(self -= self)
				.def(self * self)
				.def(self *= self)
				.def(self / self)
				.def(self /= self)


				//int
				.def(self += other<int>())
				.def(self + other<int>())
				.def(other<int>() + self)

				.def(self -= other<int>())
				.def(self - other<int>())
				.def(other<int>() - self)
				
				.def(other<int>() * self)
				.def(self * other<int>())
				.def(self *= other<int>())

				.def(other<int>() / self)
				.def(self / other<int>())
				.def(self /= other<int>())







				.def(self < self)
				.def(self <= self)
				.def(self > self)
				.def(self >= self)

				.def(self == self)

				.def(pow(self, other<bmp>() ))
				.def(pow(self, other<int>() ))

				// .def(self_ns::str(self))

				.def(abs(self))
				

				// .def("default_precision", (unsigned (bmp::*)())&bmp::default_precision) 
				// .def("default_precision", (void (bmp::*)(unsigned int))&bmp::default_precision).staticmethod("default_precision")
				.def("precision", bmpprec1)
				.def("precision", bmpprec2)

//			    .staticmethod("default_precision")

			    .def("__str__", &detail::bmp_as_str)
			    .def("__repr__", &detail::full_precision_string)
			    // .def("__format__", &detail::bmp_format)
			;
			
			def("default_precision", def_prec1);
			def("default_precision", def_prec2);


			def("sqrt",&detail::bmpsqrt);


			def("abs",&abs);
			def("sqrt",&sqrt);

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


