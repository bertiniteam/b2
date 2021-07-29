//This file is part of Bertini 2.
//
//have_bertini.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//have_bertini.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with have_bertini.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, University of Wisconsin Eau Claire


/**
\file have_bertini.hpp 

\brief Declares a function that one can use to test whether Bertini2 exists in library form.
*/


#ifndef BERTINI2_HAVE_BERTINI_HPP
#define BERTINI2_HAVE_BERTINI_HPP


/**
\brief Check for presence of the Bertini 2 library.

This function's sole purpose is for checking for the presence of the Bertini 2 library.

\return The character 'y'.
*/
extern "C"
char HaveBertini2();

#endif

