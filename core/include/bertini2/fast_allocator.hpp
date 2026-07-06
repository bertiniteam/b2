//This file is part of Bertini 2.
//
//bertini2/fast_allocator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/fast_allocator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/fast_allocator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/fast_allocator.hpp

\brief Optional faster allocator for GMP/MPFR/MPC limb allocation.

Adaptive-precision tracking spends roughly half its time in mpfr/mpc allocate/free of small,
short-lived limb buffers.  Routing GMP's allocator (which mpfr and mpc allocate through) to a
fast thread-caching allocator (mimalloc) makes each of those allocations far cheaper, with no
change to the numerics -- measured ~16% wall-clock on a representative adaptive solve.

This installs an ownership-aware hook via GMP's `mp_set_memory_functions`, calling mimalloc's
`mi_*` explicitly (mimalloc is built NON-override, so the process `malloc` is untouched -- safe
for a Python extension sharing the interpreter with NumPy and other GMP users such as gmpy2).
The hook is ownership-aware (`mi_is_in_heap_region`) so it correctly frees/reallocs GMP objects
that were allocated by the default allocator before the hook was installed (e.g. during static
init, or by another library imported first) -- so it is safe to install at any time.

When built with `-DBERTINI2_FAST_ALLOC=OFF` (no mimalloc), `InstallFastAllocator()` is a no-op.
At runtime, setting the environment variable `BERTINI2_NO_FAST_ALLOC` also disables it.
*/

#ifndef BERTINI_FAST_ALLOCATOR_HPP
#define BERTINI_FAST_ALLOCATOR_HPP

namespace bertini {

/**
\brief Install the fast multiprecision allocator (idempotent; safe to call more than once).

Call once, as early as possible, before heavy multiprecision work: at the top of `main` in the
CLI, and at module import in the Python bindings.  No-op if built without mimalloc or if the
environment variable `BERTINI2_NO_FAST_ALLOC` is set.
*/
void InstallFastAllocator();

} // namespace bertini

#endif // BERTINI_FAST_ALLOCATOR_HPP
