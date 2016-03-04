
#include "bertini_python.hpp"


namespace bertini
{
	namespace python
	{
		


		BOOST_PYTHON_MODULE(libpybertini) // this name must match the name of the generated .so file.
		{
			ExportContainers();
			
			ExportMpfr();
			
			ExportMinieigen();
		
			SetupFunctionTree();
			ExportNode();
			ExportSymbols();
			ExportOperators();
			ExportRoots();
			
			ExportSystem();
			
			ExportParsers();
		}
	
	}
}


