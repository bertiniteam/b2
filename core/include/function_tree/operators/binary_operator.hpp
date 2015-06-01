//This file is part of Bertini 2.0.
//
//Foobar is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Foobar is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// binary_operator.hpp:  Declares the class BinaryOperator.



#ifndef b2Test_BinaryOperator_h
#define b2Test_BinaryOperator_h

#include <cmath>
#include <vector>

#include "function_tree/node.hpp"

namespace bertini {
	
	// Description: This class is an interface for all binary operators, such as division.
	// Children of the operator are stored in a vector.
	class BinaryOperator : public virtual Operator
	{
	public:
		
		
		virtual void Reset() = 0; // override



	protected:
		
		virtual void print(std::ostream & target) const = 0;
				

	};

} // re: namespace bertini




#endif
