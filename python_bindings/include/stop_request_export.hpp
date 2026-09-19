//This file is part of Bertini 2.
//
//python/stop_request_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/stop_request_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/stop_request_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

//  python/stop_request_export.hpp:  exposing the cooperative stop request to Python.

#pragma once
#ifndef BERTINI_PYTHON_STOP_REQUEST_EXPORT_HPP
#define BERTINI_PYTHON_STOP_REQUEST_EXPORT_HPP

#include "python_common.hpp"


namespace bertini{
    namespace python{

        using namespace boost::python;

        /// Bind request_stop / clear_stop_request / stop_requested at module level.
        void ExportStopRequest();

}} // namespaces



#endif // include guard
