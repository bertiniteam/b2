#include "export_mpfr.hpp"
#include "mpfr_visitors.hpp"

#include <bertini2/mpfr_complex.hpp>





namespace bertini {
	namespace python{

		using namespace boost::python;


		void ExportMpfr(){
			class_<bmp>("mpfr_float", init<>())
			.def(MPFRFloatVisitor<bmp>());
			
			class_<bertini::complex>("mpfr_complex", init<>())
			.def(MPFRComplexVisitor<bertini::complex>());

		}

	}
}


