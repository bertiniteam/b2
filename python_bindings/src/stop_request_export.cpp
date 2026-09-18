//This file is part of Bertini 2.
//
//python/stop_request_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/stop_request_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/stop_request_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

//  python/stop_request_export.cpp:  the cooperative stop request, and the Ctrl-C watch that
//  drives it during a solve.

#include "stop_request_export.hpp"

#include <atomic>
#include <csignal>

namespace bertini{
    namespace python{

        namespace {

            /// Set by the handler, read by the watch after the solve returns.  Separate from the
            /// core stop flag on purpose: a stop can be requested programmatically, from another
            /// thread, and that should NOT be reported to the user as a KeyboardInterrupt.
            std::atomic<bool> g_interrupt_fired{false};
            static_assert(decltype(g_interrupt_fired)::is_always_lock_free,
                          "set from a signal handler, so it must be lock-free");

            /// The handler.  Nothing here may allocate, lock, or touch Python: two lock-free
            /// atomic stores are the whole job.
            extern "C" void OnInterrupt(int)
            {
                g_interrupt_fired.store(true, std::memory_order_release);
                bertini::RequestStop();
            }

        } // namespace

        ScopedInterruptWatch::ScopedInterruptWatch()
        {
            g_interrupt_fired.store(false, std::memory_order_release);
#ifdef _WIN32
            previous_ = std::signal(SIGINT, OnInterrupt);
#else
            // sigaction rather than signal(): on some POSIX systems signal() has one-shot
            // semantics, so a second Ctrl-C would kill the process outright instead of being
            // noticed.  No SA_RESTART: nothing in the solve is a blocking syscall we want resumed.
            struct sigaction now {};
            now.sa_handler = OnInterrupt;
            sigemptyset(&now.sa_mask);
            now.sa_flags = 0;
            sigaction(SIGINT, &now, &previous_);
#endif
        }

        ScopedInterruptWatch::~ScopedInterruptWatch()
        {
#ifdef _WIN32
            std::signal(SIGINT, previous_);
#else
            sigaction(SIGINT, &previous_, nullptr);
#endif
            // clean_ is destroyed after this body runs, withdrawing the stop request last.
        }

        bool ScopedInterruptWatch::Fired()
        {
            return g_interrupt_fired.load(std::memory_order_acquire);
        }


        void ExportStopRequest()
        {
            def("request_stop", &bertini::RequestStop,
                "Ask the solve that is running right now to stop at its next step.\n\n"
                "Cooperative and process-wide: every path being tracked notices between steps, "
                "returns SuccessCode.ExternallyTerminated, and the solve unwinds normally.  Paths "
                "that finished stay finished and inspectable; paths in flight are abandoned with "
                "where they got to; paths not yet begun are left NeverStarted and are not "
                "recorded.  Call this from another thread while solve() is blocking.  A request "
                "made while no solve is running is discarded when the next one starts -- it "
                "means 'stop what is happening', not 'do not start'.\n\n"
                "Pressing Ctrl-C during a solve does the same thing and then raises "
                "KeyboardInterrupt; this function does not raise, it just stops.");
            def("clear_stop_request", &bertini::ClearStopRequest,
                "Withdraw a stop request.  solve() does this itself on the way in and on the way "
                "out, so you only need it if you are driving a bare tracker by hand.");
            def("stop_requested", &bertini::StopRequested,
                "Has somebody asked this process to stop?  What a tracker reads between steps.");
        }

}} // namespaces
