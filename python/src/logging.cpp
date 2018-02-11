
//This file is part of Bertini 2.
//
//python/src/logging.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/src/logging.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/src/logging.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include"
//
//  Dani Brake
//  UWEC
//  Spring 2018
//
//
//  python/src/logging.cpp:  source file for exposing logging to python.


#include "logging.hpp"

namespace bertini{ namespace python{ 

void ExportSeverityLevels()
{
	enum_<logging::severity_level>("severity_level")
		.value("Debug", logging::severity_level::debug)
		.value("Trace", logging::severity_level::trace)
		.value("Info", logging::severity_level::info)
		.value("Warning", logging::severity_level::warning)
		.value("Error", logging::severity_level::error)
		.value("Fatal", logging::severity_level::fatal)
		;
}


void ExportLogging()
{
	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".logging");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("logging") = new_submodule;

	scope new_submodule_scope = new_submodule;


	ExportSeverityLevels();


	def("init", &bertini::logging::Logging::Init, 
		(boost::python::arg("pattern") = "pybertini_%N.log", boost::python::arg("format") = "%Message%", boost::python::arg("rotation_size") = 10*1024*1024, boost::python::arg("level") = logging::severity_level::info), "Initialize logging. See set_level and add_file.");

	def("set_level", &bertini::logging::Logging::SetLevel, (boost::python::arg("level")), "Set the threshold severity level.  Events with lower-than-this level will be ignored.  All messages are written to files.  Writing to strings back into Python is not currently enabled.  If this is a problem, please file an issue on GitHub at github.com/bertiniteam/b2/issues .  YAGNI");

	def("add_file", &bertini::logging::Logging::AddFile, 
		((boost::python::arg("pattern")), boost::python::arg("format"), boost::python::arg("rotation_size")), 
		"Add a file-name pattern to be written to, with a given formatting, and a threshold rotation size.  See Boost.Log for more information on these strings.  This part of PyBertini is a direct shunt to Boost.Log.");
	


}



}} // namespaces
