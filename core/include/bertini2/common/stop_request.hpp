//This file is part of Bertini 2.
//
//bertini2/common/stop_request.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/common/stop_request.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/common/stop_request.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/common/stop_request.hpp

\brief Asking a long computation to stop, from outside it.

A solve can run for hours over thousands of paths, and there was no way to tell it to stop.
In an interactive session that is worse than slow: pressing Ctrl-C sets a flag CPython raises
from its evaluation loop, and the evaluation loop is not running, because the calling thread
is down inside the solve.  The only exit was killing the process.

The mechanism here is deliberately small: one process-wide flag that anybody may set and that
long-running code reads between units of work.  Path tracking checks it between steps (see
`bertini::tracking::Tracker::TrackPath`) and reports `SuccessCode::ExternallyTerminated`, a
value that has been in the enum, unproduced, since long before this file.

Process-wide is the right scope, not per-solve.  A keyboard interrupt is a fact about the
process, not about one object; if two solves are running in threads, Ctrl-C means stop both.

Whoever sets the flag owns clearing it.  A flag left set would stop the next computation
before its first step, so the setter should use `ScopedStopRequest` or clear it by hand.  The
Python bindings install a signal handler for the duration of a solve and do exactly that.

This is cooperative.  Nothing is killed, no exception crosses a thread; the computation
notices, unwinds the way it always does, and leaves its partial results intact and
inspectable.  Work already completed stays completed, and -- since a solve records each path
as it finishes -- re-running the same ask recalls what was done and tracks only the rest.
*/

#ifndef BERTINI2_COMMON_STOP_REQUEST
#define BERTINI2_COMMON_STOP_REQUEST

#pragma once

namespace bertini {

/**
\brief Ask any running computation in this process to stop at its next check.

Safe to call from a signal handler: it only stores to an atomic.
*/
void RequestStop();

/**
\brief Withdraw a stop request, so later computations run normally.
*/
void ClearStopRequest();

/**
\brief Has somebody asked this process to stop?

Read between units of work.  The load is relaxed: this is a request, not a synchronization
point, and being one step late to notice it costs one step.
*/
bool StopRequested();

/**
\brief Clears any stale stop request on the way in, and again on the way out.

Wrap a computation in one of these when you are the one who might set the flag, so that an
abandoned request cannot leak into whatever runs next.

\code
{
    bertini::ScopedStopRequest guard;   // starts clean
    ... install a signal handler that calls RequestStop() ...
    solver.Solve();                     // stops early if the flag gets set
}                                       // flag cleared however we leave
\endcode
*/
class ScopedStopRequest
{
public:
    ScopedStopRequest() { ClearStopRequest(); }
    ~ScopedStopRequest() { ClearStopRequest(); }

    ScopedStopRequest(ScopedStopRequest const&) = delete;
    ScopedStopRequest& operator=(ScopedStopRequest const&) = delete;
};

} // namespace bertini

#endif // include guard
