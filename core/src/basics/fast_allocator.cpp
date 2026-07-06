//This file is part of Bertini 2.
//
//src/basics/fast_allocator.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/basics/fast_allocator.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/basics/fast_allocator.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/fast_allocator.hpp"

#ifdef BERTINI2_HAVE_MIMALLOC

#include <cstdlib>
#include <cstring>
#include <mutex>

#include <gmp.h>        // mp_set_memory_functions
#include <mimalloc.h>   // mi_malloc, mi_realloc, mi_free, mi_is_in_heap_region

namespace {

// GMP's allocator interface.  GMP (and therefore MPFR and MPC) routes every limb allocation
// through these.  We forward to mimalloc, but stay OWNERSHIP-AWARE so that a pointer which was
// allocated by the default allocator before this hook was installed (static init, eigenpy, an
// already-imported gmpy2, ...) is still freed/realloced correctly.  GMP hands us the old size on
// realloc/free, which makes migrating a foreign pointer exact.

void* fast_alloc(std::size_t n)
{
	return mi_malloc(n);
}

void* fast_realloc(void* p, std::size_t old_size, std::size_t new_size)
{
	if (mi_is_in_heap_region(p))
		return mi_realloc(p, new_size);

	// Foreign pointer (allocated before the hook, or by another library): migrate into mimalloc.
	void* q = mi_malloc(new_size);
	if (q && p)
		std::memcpy(q, p, old_size < new_size ? old_size : new_size);
	std::free(p);   // p came from the system allocator
	return q;
}

void fast_free(void* p, std::size_t /*size*/)
{
	if (mi_is_in_heap_region(p))
		mi_free(p);
	else
		std::free(p);   // foreign pointer; free where it was allocated
}

} // anonymous namespace

namespace bertini {

void InstallFastAllocator()
{
	static std::once_flag once;
	std::call_once(once, []{
		if (std::getenv("BERTINI2_NO_FAST_ALLOC"))
			return;
		mp_set_memory_functions(&fast_alloc, &fast_realloc, &fast_free);
	});
}

} // namespace bertini

#else // BERTINI2_HAVE_MIMALLOC not defined: built without the fast allocator

namespace bertini {

void InstallFastAllocator() { /* no-op: built with -DBERTINI2_FAST_ALLOC=OFF */ }

} // namespace bertini

#endif // BERTINI2_HAVE_MIMALLOC
