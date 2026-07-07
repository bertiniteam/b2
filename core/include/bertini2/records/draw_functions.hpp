//This file is part of Bertini 2.
//
//draw_functions.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//draw_functions.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with draw_functions.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file draw_functions.hpp

\brief Dependency-free declarations of the pinned draw primitives (ADR-0044).

The foundational headers (num_traits, double_extensions) need double draws from the
per-thread pinned stream but must not include the full derive.hpp (which needs the
multiprecision types those headers define).  These free functions break the cycle;
definitions live in src/records/derive.cpp.
*/

#pragma once

namespace bertini {
namespace records {

/// \brief A draw from this thread's pinned stream, uniform on [0,1) (exactly 53 bits).
double DrawUnitDouble();

/// \brief A draw from this thread's pinned stream, uniform on [-1,1].
double DrawSymmetricDouble();

} // namespace records
} // namespace bertini
