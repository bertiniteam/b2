//This file is part of Bertini 2.
//
//bertini2/flint_version.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/flint_version.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/flint_version.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file flint_version.hpp

\brief The version of FLINT the library was linked against.

Part of the FLINT spike (bertiniteam/b2#412): a reference to one FLINT symbol from inside
libbertini2, so that every wheel really carries FLINT and its size can be measured.  A linker
drops a library nothing references, and delvewheel then has nothing to bundle.
*/

#pragma once

#include <string>

namespace bertini {

/**
\brief The FLINT version string the library was linked against, e.g. "3.5.0".
*/
std::string FlintVersion();

} // namespace bertini
