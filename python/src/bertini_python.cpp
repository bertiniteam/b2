//This file is part of Bertini 2.0.
//
// python/bertini_python.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/bertini_python.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/bertini_python.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/bertini_python.cpp:  the main source file for the python interface for bertini.


#include "bertini_python.hpp"


namespace bertini
{
	namespace python
	{
		


		BOOST_PYTHON_MODULE(pybertini) // this name must match the name of the generated .so file.
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


