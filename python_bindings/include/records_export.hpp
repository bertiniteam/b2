//This file is part of Bertini 2.
//
//python_bindings/include/records_export.hpp is free software: you can redistribute it
//and/or modify it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
//See the GNU General Public License for more details.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.
// Bertini2 is provided with permitted additional terms in the b2/licenses/ directory.

/** \file records_export.hpp
\brief Exports the `records` submodule: the structured output directory (b2rec/1).
*/

#pragma once
#include "python_common.hpp"

namespace bertini { namespace python {

/// \brief Export the `records` submodule (OutputDirectory and friends).
void ExportRecords();

}} // namespaces
