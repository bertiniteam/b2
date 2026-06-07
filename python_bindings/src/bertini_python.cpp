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
// Copyright(C) Bertini2 Development Team
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
//  silviana amethyst
//  UWEC
//  Spring 2018, Summer 2023
//
//
//  python/bertini_python.cpp:  the main source file for the python interface for bertini.

// Only the types/functions directly used in this TU are included.
// The mega-include bertini_python.hpp is intentionally NOT included here:
// it pulls in the full tracker/endgame/NID header stack (causing ~90s compile
// time) even though this file only calls void export functions by name.
// gather_variables registration was moved into SetupFunctionTree() in
// function_tree_export.hpp where those types are already available.

// Include only what this TU directly needs.
// The heavy tracker/endgame/NID headers are intentionally excluded:
// the export functions they define are called here by name only (void, no args).
#include "function_tree_export.hpp"  // SetupFunctionTree() + light function_tree headers
#include "eigenpy_interaction.hpp"   // EnableEigenPy()
#include "parallel_export.hpp"       // ExportParallel()

namespace bertini { namespace python {

// Forward-declare the heavy export functions — their headers pull in
// the full tracker/endgame/NID stack which is not needed in this TU.
void ExportContainers();
void ExportDetails();
void ExportMpfr();
void ExportRandom();
void ExportAllSystems();
void ExportParsers();
void ExportTrackers();
void ExportTrackerObservers();
void ExportEndgames();
void ExportEndgameObservers();
void ExportLogging();
void ExportZeroDim();
void ExportNID();
void ExportInfo();

} } // namespace bertini::python


namespace bertini
{
	namespace python
	{

		BOOST_PYTHON_MODULE(_pybertini) // this name must match the name of the generated .so file.
		{
			// see https://stackoverflow.com/questions/6114462/how-to-override-the-automatically-created-docstring-data-for-boostpython
			docstring_options docopt;
			docopt.enable_all();
			docopt.disable_cpp_signatures();

			object package = scope();
		    package.attr("__path__") = "_pybertini";

		    // do this one first, so that the later calls into EigenPy work :)
		    EnableEigenPy();

			ExportContainers();

			ExportDetails();

			ExportMpfr();

			ExportRandom();

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

				boost::python::def("gather_variables",
					static_cast<bertini::VariableGroup(*)(std::vector<std::shared_ptr<bertini::node::Function>> const&)>(&bertini::node::GatherVariables),
					(boost::python::arg("functions")),
					"Return the distinct variables appearing in a list of functions, ordered alphabetically by name.");
			}

			ExportAllSystems();

			ExportParsers();

			ExportTrackers();
			ExportTrackerObservers();

			ExportEndgames();
			ExportEndgameObservers();

			ExportLogging();

			ExportParallel();
			ExportZeroDim();

			ExportNID();

			ExportInfo();
		}

	}
}
