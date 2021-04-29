//This file is part of Bertini 2.
//
//bertini2/common/stream_enum.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/common/stream_enum.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with tracking/include/bertini2/trackers/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


// individual authors of this file include:
// silviana amethyst, university of notre dame, university of wisconsin eau claire
// Tim Hodges, Colorado State University

#ifndef BERTINI2_STREAM_ENUM
#define BERTINI2_STREAM_ENUM

#pragma once


namespace bertini
{
	
	/**
	\brief Method for printing class enums to streams.

	This was adapted from https://stackoverflow.com/questions/11421432/how-can-i-output-the-value-of-an-enum-class-in-c11
	asked by user Adi, answered by James Adkison.  This code was provided CC-BY-SA 3.

	This code does NOT work for streaming enum classes to Boost.Log streams.

	\param stream The stream to print to.
	\param e The enum value to print
	\return The stream you are writing to.
	*/
	template<typename T>
	std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e)
	{
	    return stream << static_cast<typename std::underlying_type<T>::type>(e);
	}

} // namespace bertini




#endif // include guard

