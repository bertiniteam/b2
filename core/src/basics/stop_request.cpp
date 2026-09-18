//This file is part of Bertini 2.
//
//src/basics/stop_request.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/basics/stop_request.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/basics/stop_request.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/common/stop_request.hpp"

#include <atomic>

namespace bertini {

namespace {

/// The one flag.  `sig_atomic_t`-style lock-free storage is what makes RequestStop safe to
/// call from a signal handler; a std::atomic<bool> is lock-free on every platform we build
/// for, and the static_assert below says so out loud rather than hoping.
std::atomic<bool> g_stop_requested{false};
static_assert(decltype(g_stop_requested)::is_always_lock_free,
              "the stop flag is set from a signal handler, so it must be lock-free");

} // namespace

void RequestStop()
{
    // release: everything this thread did before asking is visible to whoever sees the request
    g_stop_requested.store(true, std::memory_order_release);
}

void ClearStopRequest()
{
    g_stop_requested.store(false, std::memory_order_release);
}

bool StopRequested()
{
    // relaxed: read on a hot path, between steps of every path being tracked.  Noticing one
    // step late costs one step, and there is nothing to synchronize with -- the flag carries
    // no payload, it is the whole message.
    return g_stop_requested.load(std::memory_order_relaxed);
}

} // namespace bertini
