//This file is part of Bertini 2.0.
//
//nary_operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nary_operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nary_operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// nary_operator.hpp:  Declares the class NaryOperator.

#ifndef b2Test_NaryOperator_h
#define b2Test_NaryOperator_h

#include <vector>

#include "function_tree/node.hpp"
#include "function_tree/operators/operator.hpp"


namespace bertini {
	
	
	// Description: This class is an interface for all n-ary operators, such as summation and multiplication.
	// Children of the operator are stored in a vector and methods to add and access children are available
	// in this interface.
	class NaryOperator : public virtual Node, public virtual Operator
	{
	public:
		
		virtual ~NaryOperator() = default;
		
		
		virtual void Reset() override
		{
			Node::ResetStoredValues();
			for (auto ii:children_)
			{
				ii->Reset();
			}
		}
		
		// Add a child onto the container for this operator
		virtual void AddChild(std::shared_ptr<Node> child) // override
		{
			children_.push_back(std::move(child));
		}
		
		
		
		
		
		
		size_t children_size()
		{
			return children_.size();
		}
		
		std::shared_ptr<Node> first_child()
		{
			return children_[0];
		}
		
		
		
		////////////// TESTING /////////////////
		virtual void PrintTree() override
		{
			for(int ii = 0; ii < tabcount; ++ii)
			{
				std::cout << "\t";
			}
			std::cout << tabcount+1 << "." <<  boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>() << std::endl;
			tabcount++;
			for(auto vv : children_)
			{
				vv->PrintTree();
			}
			tabcount--;
		}
		////////////// TESTING /////////////////
		
		
		
		
		
	protected:
		//Stores all children for this operator node.
		//This is an NaryOperator and can have any number of children.
		std::vector< std::shared_ptr<Node> > children_;
		
		
		
		
	};
	
	
} // re: namespace bertini


#endif
