//This file is part of Bertini 2.
//
//python/numerical_irreducible_decomposition_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/numerical_irreducible_decomposition_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/numerical_irreducible_decomposition_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  University of Wisconsin-Eau Claire
//  2026
//
//
//  python/numerical_irreducible_decomposition_export.cpp:  Coordinator — sets up the
//  nag_algorithms submodule and calls per-tracker-family registration functions.
//  No heavy endgame/NID headers needed here: only void functions are called by name.

namespace bertini{
	namespace python{

		// Forward declarations — defined in nid_{datatypes,double,mp,amp}_export.cpp.
		void ExportNIDDataTypes();
		void ExportNIDDouble();
		void ExportNIDMP();
		void ExportNIDAMP();

		void ExportNID(){
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".nag_algorithms");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("nag_algorithms") = new_submodule;

			scope new_submodule_scope = new_submodule;

			// config structs (Tolerances, Regeneration, Sharpening, PostProcessing)
			// are registered by ExportZeroDim which runs before this.
			ExportNIDDataTypes();
			ExportNIDDouble();
			ExportNIDMP();
			ExportNIDAMP();
		}

	}
}
