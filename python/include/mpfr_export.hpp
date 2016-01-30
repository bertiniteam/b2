




#ifndef BERTINI_PYTHON_BOOST_MP_HPP
#define BERTINI_PYTHON_BOOST_MP_HPP
#include "mpfr_visitors.hpp"



namespace bertini{
	namespace python{
		
		void ExportMpfr()
		{
			class_<bmp>("mpfr_float", init<>())
			.def(MPFRFloatVisitor<bmp>());
			
			class_<bertini::complex>("mpfr_complex", init<>())
			.def(MPFRComplexVisitor<bertini::complex>());
		};

		
	} //namespace python
} // namespace bertini

#endif

