//This file is part of Bertini 2.
//
//python/python_common.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/python_common.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/python_common.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

//  python/python_common.hpp:  A common header file for all python exposure files

#pragma once
#ifndef BERTINI_PYTHON_COMMON_HPP
#define BERTINI_PYTHON_COMMON_HPP

#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python/wrapper.hpp>

#include <boost/python/operators.hpp>
#include <boost/operators.hpp>

#include <sstream>


#include <bertini2/mpfr_complex.hpp>
#include <bertini2/mpfr_extensions.hpp>
#include <bertini2/eigen_extensions.hpp>
#include <bertini2/common/stop_request.hpp>

#include <csignal>
#ifndef _WIN32
#include <signal.h>      // struct sigaction, which <csignal> alone need not expose
#endif



using namespace boost::python;
using real_mp = bertini::real_mp;
using complex_mp = bertini::complex_mp;

namespace bertini { namespace python {

/**
\brief RAII: release the Python GIL for the duration of a long C++ call (a threaded solve), so
worker threads run truly in parallel and Python observer callbacks can re-acquire the GIL.
Equivalent to Py_BEGIN_ALLOW_THREADS / Py_END_ALLOW_THREADS.
*/
struct ScopedGILRelease
{
    PyThreadState* state_;
    ScopedGILRelease()  : state_(PyEval_SaveThread()) {}
    ~ScopedGILRelease() { PyEval_RestoreThread(state_); }
    ScopedGILRelease(ScopedGILRelease const&) = delete;
    ScopedGILRelease& operator=(ScopedGILRelease const&) = delete;
};

/**
\brief RAII: ensure the calling thread holds the GIL before touching Python objects, and restore
on scope exit.  Safe to call from any thread -- a C++ worker thread (the solve released the GIL)
or the main thread.  Used at the top of the observer trampoline before calling into a Python
observer's Observe().
*/
struct ScopedGILAcquire
{
    PyGILState_STATE state_;
    ScopedGILAcquire()  : state_(PyGILState_Ensure()) {}
    ~ScopedGILAcquire() { PyGILState_Release(state_); }
    ScopedGILAcquire(ScopedGILAcquire const&) = delete;
    ScopedGILAcquire& operator=(ScopedGILAcquire const&) = delete;
};

/**
\brief RAII: for the duration of a long C++ call, make Ctrl-C stop it instead of being ignored.

Why this exists.  A solve releases the GIL and runs on the calling thread.  When the user
presses Ctrl-C, CPython's C-level handler notes the signal and waits for the evaluation loop
to raise KeyboardInterrupt -- and the evaluation loop is not running, because the thread is
inside the solve.  So the interrupt is delivered, recorded, and ignored until the very thing
the user wanted to stop has finished.  Killing the process was the only exit.

What this does.  For the scope of the object it installs its own SIGINT handler, which does
two async-signal-safe things: calls bertini::RequestStop(), which every tracking loop checks
between steps, and marks that the signal fired.  The solve then unwinds the way it always
does, cooperatively, leaving completed paths installed and inspectable.  On scope exit the
previous handler is restored -- Python's, or whatever the user had -- and the stop request is
withdrawn so nothing leaks into the next solve.  The caller reads Fired() and, if so, raises
KeyboardInterrupt in Python, which is what the user asked for.

Only the C-level handler is touched.  Python's record of the Python-level handler
(signal.getsignal) is never changed, so a user's own handler comes back exactly as it was.

Not applied around the MPI path: a signal reaches one rank, and stopping one rank of a
manager-worker solve cleanly is its own problem.
*/
class ScopedInterruptWatch
{
public:
    ScopedInterruptWatch();
    ~ScopedInterruptWatch();
    ScopedInterruptWatch(ScopedInterruptWatch const&) = delete;
    ScopedInterruptWatch& operator=(ScopedInterruptWatch const&) = delete;

    /// Did SIGINT arrive while this watch was installed?
    static bool Fired();

private:
    bertini::ScopedStopRequest clean_;   ///< declared first, so destroyed last: the flag is withdrawn after the handler is gone
#ifdef _WIN32
    void (*previous_)(int) = nullptr;
#else
    struct sigaction previous_ {};
#endif
};

}} // namespace bertini::python

#endif
