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
// unary_operator.hpp:  Declares the class UnaryOperator.



#ifndef b2Test_UnaryOperator_h
#define b2Test_UnaryOperator_h


#include <vector>

#include "function_tree/node.hpp"


namespace bertini {
	
	// Description: This class is an interface for all unary operators, such as negation.
	// The sole child is stored in a shared_ptr.
	class UnaryOperator : public virtual Node, public virtual Operator
	{
	public:
		
		UnaryOperator(){}
		
		UnaryOperator(const std::shared_ptr<Node> & N) : child_(N)
		{}
		
		
		
		virtual ~UnaryOperator() = default;
		
		
		virtual void Reset() override
		{
			Node::ResetStoredValues();
			child_->Reset();
		}
		
		
		
		virtual void SetChild(std::shared_ptr<Node> new_child)
		{
			child_ = new_child;
		}
		
		
		
		//Return the only child for the unary operator
		std::shared_ptr<Node> first_child()
		{
			return child_;
		}
		
		
		////////////// TESTING /////////////////
		virtual void PrintTree()
		{
			for(int ii = 0; ii < tabcount; ++ii)
			{
				std::cout << "\t";
			}
			std::cout << tabcount+1 << "." <<  boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
			tabcount++;
			child_->PrintTree();
			tabcount--;
		}
		////////////// TESTING /////////////////
		
		
		
		
		
		
	protected:
		//Stores the single child of the unary operator
		std::shared_ptr<Node> child_;
	};
	
} // re: namespace bertini


#endif
