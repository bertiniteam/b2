//This file is part of Bertini 2.
//
//python/info.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/info.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/info.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

//  python/info.cpp:  source file for exposing trackers to python.


#include "info.hpp"

#include "bertini2/flint_version.hpp"

namespace bertini{
    namespace python{

        void ExportInfo()
        {
            scope current_scope;
            std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
            new_submodule_name.append(".info");
            object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
            current_scope.attr("info") = new_submodule;

            scope new_submodule_scope = new_submodule;
            new_submodule_scope.attr("__doc__") = "Meta information about the library";


            def("version", &bertini::Version);
            def("boost_version", &bertini::BoostHeaderVersion);
            def("eigen_version", &bertini::EigenHeaderVersion);
            def("gmp_version", &bertini::GMPVersion);
            def("mpfr_version", &bertini::MPFRVersion);
            // Also the reference that makes the extension depend on FLINT: libbertini2 is a static
            // library, so the linker drops flint_version.o, and with it the FLINT import the
            // wheel repair tools bundle, unless something here uses it.
            def("flint_version", &bertini::FlintVersion);
            def("mpi_version", &bertini::MPIVersion);
            def("dependency_versions", &bertini::DependencyVersions);
            def("splash", &bertini::SplashScreen);
        }
    }
}
