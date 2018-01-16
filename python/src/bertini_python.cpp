//This file is part of Bertini 2.
//
//python/bertini_python.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/bertini_python.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/bertini_python.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//  Danielle Brake
//  UWEC
//  Spring 2018
//
//
//  python/bertini_python.cpp:  the main source file for the python interface for bertini.


#include "bertini_python.hpp"


namespace bertini
{
	namespace python
	{
		


		BOOST_PYTHON_MODULE(_pybertini) // this name must match the name of the generated .so file.
		{
			// see https://stackoverflow.com/questions/6114462/how-to-override-the-automatically-created-docstring-data-for-boostpython
			// docstring_options d(true, true, false); // local_
			docstring_options docopt;
			docopt.enable_all();
			docopt.disable_cpp_signatures();

			object package = scope();
		    package.attr("__path__") = "_pybertini";

			ExportContainers();
			
			ExportDetails();

			ExportMpfr();
			
			ExportMinieigen();
		
			SetupFunctionTree();

			{
				scope current_scope;
				

				std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
				new_submodule_name.append(".function_tree");
				object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
				current_scope.attr("function_tree") = new_submodule;
				

				scope new_submodule_scope = new_submodule;
				new_submodule_scope.attr("__doc__") = "The symbolics for Bertini2.  Operator overloads let you write arithmetic do form your system, after making variables, etc.";
				ExportNode();
				ExportSymbols();
				ExportOperators();
				ExportRoots();
			}

			ExportAllSystems();
			
			ExportParsers();

			ExportTrackers();
			ExportTrackerObservers();

			ExportEndgames();
			ExportEndgameObservers();

			
		}
	
	}
}


