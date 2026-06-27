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

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
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



using namespace boost::python;
using mpfr_float = bertini::mpfr_float;
using mpfr_complex = bertini::mpfr_complex;

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

}} // namespace bertini::python

#endif
