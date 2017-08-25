//This file is part of Bertini 2.
//
//bertini2/parallel/initialize_finalize.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/initialize_finalize.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/initialize_finalize.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/parallel/initialize_finalize.hpp 

\brief Provides the free initialization and finalization routines for Bertini2.
*/


#pragma once

#include <iostream>

#include "bertini2/io/splash.hpp"

namespace bertini{

	namespace parallel{

	void Finalize();
	void Initialize();	

	}

	namespace serial{
		void Finalize();
		void Initialize();
	}
}
