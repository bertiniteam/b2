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

#ifndef operator_h
#define operator_h

#include "function_tree/node.hpp"

#include <vector>


namespace bertini {
	
	// Description: This class is could work as an interface for all operators in a function tree.
	class Operator : public virtual Node
	{
		
	public:
		
		virtual ~Operator() = default;
	};
	
	
	
	
	// Description: This class is an interface for all unary operators, such as negation.
	// The sole child is stored in a shared_ptr.
	class UnaryOperator : public virtual Operator
	{
	public:
		
		UnaryOperator(){}
		
		UnaryOperator(const std::shared_ptr<Node> & N) : child_(N)
		{}
		
		
		
		virtual ~UnaryOperator() = default;
		
		
		void Reset() override
		{
			Node::ResetStoredValues();
			child_->Reset();
		}
		
		
		
		void SetChild(std::shared_ptr<Node> new_child)
		{
			child_ = new_child;
		}
		
		
		
		//Return the only child for the unary operator
		std::shared_ptr<Node> first_child()
		{
			return child_;
		}
		
		
		
		/**
		 Compute the degree of a node.  For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return nan.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return child_->Degree(v);
		}
		

		int Degree(std::vector< std::shared_ptr<Variable > > const& vars) const override
		{
			auto multideg = MultiDegree(vars);
			auto deg = 0;
			std::for_each(multideg.begin(),multideg.end(),[&](int n){
							if (n < 0)
								deg = n;
							else
	                        	deg += n;
	 						});
			return deg;
		}


		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	    */
		std::vector<int> MultiDegree(std::vector< std::shared_ptr<Variable> > const& vars) const override
		{
			
			std::vector<int> deg(vars.size());
			for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
			{
				*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
			}
			return deg;
		}


		void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			child_->Homogenize(vars, homvar);
		}
		
		bool IsHomogeneous() const override
		{
			if (Degree()==0)
			{
				return true;
			}
			else
				return false;
		}
		
		
	protected:
		//Stores the single child of the unary operator
		std::shared_ptr<Node> child_;
	};
	
	
	
	// Description: This class is an interface for all binary operators, such as division.
	// Children of the operator are stored in a vector.
	class BinaryOperator : public virtual Operator
	{
	public:
		
		
		virtual void Reset() = 0; // override
		
		
		
	protected:
		
		virtual void print(std::ostream & target) const = 0;
		
		
	};
	
	
	// Description: This class is an interface for all n-ary operators, such as summation and multiplication.
	// Children of the operator are stored in a vector and methods to add and access children are available
	// in this interface.
	class NaryOperator : public virtual Node, public virtual Operator
	{
	public:
		
		virtual ~NaryOperator() = default;
		
		
		void Reset() override
		{
			Node::ResetStoredValues();
			for (auto ii:children_)
			{
				ii->Reset();
			}
		}
		
		// Add a child onto the container for this operator
		virtual void AddChild(std::shared_ptr<Node> child)
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
		
		
		
		
		
		
		
	protected:
		//Stores all children for this operator node.
		//This is an NaryOperator and can have any number of children.
		std::vector< std::shared_ptr<Node> > children_;
		
		
		
		
	};
	
	
	
} // re: namespace bertini


#endif
