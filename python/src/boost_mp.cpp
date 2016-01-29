#include "boost_mp.hpp"
#include "mpfr_visitors.hpp"





namespace bertini {
	namespace python{

		using namespace boost::python;


		void ExportMpfrFloat(){
			class_<bmp>("mpfr_float", init<>())
			.def(MPFRFloatVisitor<bmp>());
		}

	}
}


