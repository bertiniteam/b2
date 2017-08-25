//This file is part of Bertini 2.
//
//limbo.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//limbo.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with limbo.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, University of Wisconsin Eau Claire


/**
\file limbo.hpp 

\brief Declares functions which have not yet found a home. 

If you can find a better home for these functions, please make the changes and submit a pull request.
*/


#ifndef BERTINI_LIMBO_HPP
#define BERTINI_LIMBO_HPP

#include <vector>
#include <stdexcept>

extern "C" {
	/**
	\brief Check for presence of the Bertini 2 library.

	This function's sole purpose is for checking for the presence of the Bertini 2 library.

	\return The character 'y'.
	*/
	char HaveBertini2();
}









namespace bertini{

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

	/**
	\brief Convert a zero-based index to a zero-based subscript vector.  

	Throws `std::out_of_range` if the index is out-of-range based on the dimensions.

	This goes from front to back, top to bottom.  So 

	[0 2 4      [(0,0) (0,1) (0,2)
	 1 3 5]      (1,0) (1,1) (1,2)]

	 etc for higher-dimensional arrays.  The functionality here is identical to that of Matlab's analagous call.
	
	\param index The index you want to convert.
	\param dimensions The dimensions of the object you are subscripting or indexing into.
	\return A vector containing the subscripts for the input index.
	\throws std::out_of_range, if the index is impossible to convert.
	*/
	template <typename T>
	std::vector<T> IndexToSubscript(T index, std::vector<T> const& dimensions)
	{

		std::vector< T > subscripts(dimensions.size());//for forming a subscript from an index

		std::vector<T> k(dimensions.size(),1);
		for (int ii = 0; ii < dimensions.size()-1; ++ii)
		  k[ii+1] = k[ii]*dimensions[ii];


		if (index >= k.back()*dimensions.back())
		  throw std::out_of_range("in IndexToSubscript, index exceeds max based on dimension sizes");


		for (int ii = dimensions.size()-1; ii >= 0; --ii)
		{
		  T I = index%k[ii];
		  T J = (index - I) / k[ii];
		  subscripts[ii] = J;
		  index = I;
		}

		return subscripts;
	  }

	  
	  /**
	  \brief Three-argument form of `max`.

	  \param a Input one
	  \param b Input two
	  \param c Input three
	  */
	  template <typename T>
	  T max(T const& a, T const& b, T const& c)
	  {
	  	using std::max;
	  	return max(max(a,b),c);
	  }

	  /**
	  \brief Four-argument form of `max`.

	  \param a Input one
	  \param b Input two
	  \param c Input three
	  \param d Input four
	  */
	  template <typename T>
	  T max(T const& a, T const& b, T const& c, T const& d)
	  {
	  	using std::max;
	  	using bertini::max;
	  	return max(max(a,b,c),d);
	  }

	  /**
	  \brief Five-argument form of `max`.

	  \param a Input one
	  \param b Input two
	  \param c Input three
	  \param d Input four
	  \param e Input five
	  */
	  template <typename T>
	  T max(T const& a, T const& b, T const& c, T const& d, T const& e)
	  {
	  	using std::max;
	  	using bertini::max;
	  	return max(max(a,b,c,d),e);
	  }
} // re: namespace bertini


#endif

