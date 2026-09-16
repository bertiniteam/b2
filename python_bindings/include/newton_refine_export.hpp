//This file is part of Bertini 2.
//
//python/newton_refine_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/newton_refine_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/newton_refine_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// silviana amethyst, university of wisconsin eau claire

/**
\file python/newton_refine_export.hpp

\brief Exposes the standalone NewtonRefine (tracker-free Newton against a System,
deflated systems above all) to Python.
*/

#pragma once

#include "python_common.hpp"

namespace bertini {
namespace python {

/**
\brief Export the standalone newton_refine free function to the current module scope.
*/
void ExportNewtonRefine();

} // namespace python
} // namespace bertini
