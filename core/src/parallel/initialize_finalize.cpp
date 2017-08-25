//This file is part of Bertini 2.
//
//bertini2/parallel/initialize_finalize.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/initialize_finalize.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/initialize_finalize.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/parallel/initialize_finalize.cpp 

\brief Provides the free initialization and finalization routines for Bertini2.
*/



#include "bertini2/parallel/initialize_finalize.hpp"


namespace bertini{

	namespace parallel{
		void Initialize()
		{

		}

		void Finalize()
		{

		}
	}

	namespace serial{
		void Initialize()
		{
			std::cout << "\n\n\n" << SplashScreen() << "\n\n\n";
			std::cout << "\n\n" << DependencyVersions() << "\n\n";
		}

		void Finalize()
		{

		}
	}

}






