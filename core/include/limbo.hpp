//This file is part of Bertini 2.0.
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
//  Created by Daniel Brake, July 13, 2015.
//
//
// limbo.hpp:  Declares functions which have not yet found a home.


#include <vector>



namespace bertini{

	/**
	Convert a zero-based index to a zero-based subscript vector.  Throws if the index is out-of-range based on the dimensions.

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

}