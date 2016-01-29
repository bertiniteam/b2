
#include "bertini_python.hpp"

BOOST_PYTHON_MODULE(libpybertini) // this name must match the name of the generated .so file.
{

	using namespace bertini::python;
	ExportMpfr();
//	ExportMpfrComplex();
	
	SetupFunctionTree();
	ExportNode();
	ExportSymbols();
	ExportOperators();
	ExportRoots();
	ExportSystem();
	
}


