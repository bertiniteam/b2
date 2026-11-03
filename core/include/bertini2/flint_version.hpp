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

FLINT 3 carries Arb, the ball arithmetic behind b2's interval types (bertiniteam/b2#412).  This
is also the one FLINT reference the Python extension makes, which keeps FLINT linked into it:
libbertini2 is a static library, so a linker drops what nothing references, and the wheel repair
tools then bundle nothing.
*/

#pragma once

#include <string>

namespace bertini {

/**
\brief The FLINT version string the library was linked against, e.g. "3.5.0".
*/
std::string FlintVersion();

} // namespace bertini
