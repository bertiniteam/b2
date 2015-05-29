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
// sum_operator.h:  Declares the class SumOperator.

#ifndef b2Test_SumOperator_h
#define b2Test_SumOperator_h

#include <vector>
#include <string>   

#include "function_tree/node.h"
#include "function_tree/operators/nary_operator.h"

namespace bertini {

// Node -> NaryOperator -> SumOperator
// Description: This class represents summation and difference operators.  All children are terms and are stored
// in a single vector, and a vector of bools is used to determine the sign of each term.  FreshEval method
// is defined for summation and difference.
class SumOperator : public virtual NaryOperator
{
public:
	virtual ~SumOperator() = default;
	
	
	
	SumOperator(){}
	
	SumOperator(const std::shared_ptr<Node> & left, const std::shared_ptr<Node> & right)
	{
		AddChild(left);
		AddChild(right);
	}
	
	
	SumOperator(const std::shared_ptr<Node> & left, bool add_or_sub_left, const std::shared_ptr<Node> & right, bool add_or_sub_right)
	{
		AddChild(left, add_or_sub_left);
		AddChild(right, add_or_sub_right);
	}
	
	
	SumOperator& operator+=(const std::shared_ptr<Node> & rhs)
	{
		this->AddChild(rhs);
		return *this;
	}
	
	SumOperator& operator-=(const std::shared_ptr<Node> & rhs)
	{
		this->AddChild(rhs,false);
		return *this;
	}
	
	
    // Print the data for this node, as well as all it's children
    //TODO (JBC): Implement this method
    virtual std::string PrintNode() override {return "";}
    
    
    
    //Special Behaviour: by default all terms added are positive
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        NaryOperator::AddChild(std::move(child));
        children_sign_.push_back(true);
    }
    
    
    
    //Special Behaviour: Pass bool to set sign of term: true = add, false = subtract
	void AddChild(std::shared_ptr<Node> child, bool sign) // not an override
    {
        NaryOperator::AddChild(std::move(child));
        children_sign_.push_back(sign);
    }


	virtual void print(std::ostream & target) const override
	{
		target << "(";
		for (auto iter = children_.begin(); iter!= children_.end(); iter++) {
			if (iter==children_.begin()) {
				// on the first iteration, no need to put a + if a +
				if ( !(*(children_sign_.begin()+(iter-children_.begin()))) )
					target << "-";
			}
			else
			{
				if ( !(*(children_sign_.begin()+(iter-children_.begin()))) )
					target << "-";
				else
					target << "+";
			}
			(*iter)->print(target);

		}
		target << ")";
	}

protected:
    // Specific implementation of FreshEval for add and subtract.
    //  If child_sign_ = true, then add, else subtract
    virtual dbl FreshEval(dbl) override
    {
        dbl retval{0};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_sign_[ii])
            {
                retval += children_[ii]->Eval<dbl>();
            }
            else
            {
                retval -= children_[ii]->Eval<dbl>();
            }
        }
        
        return retval;
    }
    
    
    
    
    virtual mpfr FreshEval(mpfr) override
    {
        mpfr retval{0};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_sign_[ii])
            {
                retval += children_[ii]->Eval<mpfr>();
            }
            else
            {
                retval -= children_[ii]->Eval<mpfr>();
            }
        }
        
        return retval;
    }

    
    
    
    
    
private:
    // Stores the sign of the particular term.  There is a one-one
    // correspondence between elements of children_sign_ and children_.  This
    // is enforced by the AddChild method below, redefined in SumOperator.
    
    // TODO(JBC): If we add method to delete child, must also delete children_sign_ entry.
    std::vector<bool> children_sign_;

};



} // re: namespace bertini


namespace  {

using Node = bertini::Node;
using SumOperator = bertini::SumOperator;

inline std::shared_ptr<Node>& operator+=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
{
	std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
	
	std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
	std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(rhs);
	
	lhs.swap(temp);
	return lhs;
}


inline std::shared_ptr<Node>& operator-=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
{
	std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
	
	std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
	std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(rhs,false);
	
	lhs.swap(temp);
	return lhs;
}


inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
{
	return std::make_shared<SumOperator>(lhs,rhs);
}


inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
{
	return std::make_shared<SumOperator>(lhs,true,rhs,false);
}


}
#endif
