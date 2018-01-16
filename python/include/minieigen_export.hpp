//This file is part of Bertini 2.
//
//python/include/minieigen_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/include/minieigen_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/include/minieigen_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
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
//  Danielle Brake
//  University of Wisconsin - Eau Claire
//  Spring 2018
//
//
//
//  python/include/minieigen_export.hpp:  Header file for exposing the eigen data types needed for bertini.



#pragma once
#ifndef BERTINI_PYTHON_MINIEIGEN_EXPORT_HPP
#define BERTINI_PYTHON_MINIEIGEN_EXPORT_HPP

#include <bertini2/function_tree.hpp>
#include <bertini2/eigen_extensions.hpp>

#include "python_common.hpp"

#include "minieigen/src/common.hpp"
#include "minieigen/src/converters.hpp"
#include "minieigen/src/visitors.hpp"





namespace bertini{
	namespace python{
		
		
		void ExportMinieigen();
		
		
	}
}


#endif
