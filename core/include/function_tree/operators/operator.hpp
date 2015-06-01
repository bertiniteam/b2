//This file is part of Bertini 2.0.
//
//operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// operator.hpp:  Declares the class Operator.

#ifndef b2Test_Operator_h
#define b2Test_Operator_h

#include "function_tree/node.hpp"


namespace bertini {
	
	// Description: This class is could work as an interface for all operators in a function tree.
	class Operator : public virtual Node
	{
		
	public:
		
		// This method adds a child to an operator.
		//
		// Every interface for operator types(nary, binary and unary) implement this method.
		//	virtual void AddChild(std::shared_ptr<Node> child) = 0;
		//
		virtual ~Operator() = default;
	};
	
} // re: namespace bertini


#endif
