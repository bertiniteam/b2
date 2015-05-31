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
// mult_operator.hpp:  Declares the class MultOperator.



#ifndef b2Test_MultOperator_h
#define b2Test_MultOperator_h

#include "function_tree/node.hpp"
#include "function_tree/operators/nary_operator.hpp"


namespace bertini {

// Node -> NaryOperator -> MultOperator
// Description: This class represents multiplication operator.  All children are factors and are stored
// in a vector.  FreshEval method is defined for multiplication.
class MultOperator : public virtual NaryOperator
{
public:
//
	
	MultOperator(){}
	
	MultOperator(const std::shared_ptr<Node> & left, const std::shared_ptr<Node> & right)
	{
		AddChild(left);
		AddChild(right);
	}
	
	
	MultOperator(const std::shared_ptr<Node> & left, bool mult_or_div_left, const std::shared_ptr<Node> & right, bool mult_or_div_right)
	{
		AddChild(left, mult_or_div_left);
		AddChild(right, mult_or_div_right);
	}
	
	
	
	
	
	
	
	
	virtual ~MultOperator() = default;
	
    // See node.hpp for description
    // TODO:(JBC) Implement this for multiplication
    virtual std::string PrintNode() override {return "";}


    //Special Behaviour: by default all factors are in numerator
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        NaryOperator::AddChild(std::move(child));
        children_mult_or_div_.push_back(true);
    }
    
    
    
    //Special Behaviour: Pass bool to set sign of term: true = mult, false = divide
	void AddChild(std::shared_ptr<Node> child, bool mult) // not an override
    {
        NaryOperator::AddChild(std::move(child));
        children_mult_or_div_.push_back(mult);
    }


	/**
	 overridden method for printing to an output stream
	 */
	virtual void print(std::ostream & target) const override
	{
		target << "(";
		for (auto iter = children_.begin(); iter!= children_.end(); iter++) {
			(*iter)->print(target);
			if (iter!=(children_.end()-1)){
				if (*(children_mult_or_div_.begin() + (iter-children_.begin())+1)) {
					target << "*";
				}
				else{
					target << "/";
				}
				
			}
		}
		target << ")";
	}
		


protected:
    // Specific implementation of FreshEval for mult and divide.
    //  If child_mult_ = true, then multiply, else divide
    virtual dbl FreshEval(dbl) override
    {
        dbl retval{1};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_mult_or_div_[ii])
            {
                retval *= children_[ii]->Eval<dbl>();
            }
            else
            {
                retval /= children_[ii]->Eval<dbl>();
            }
        }
        
        return retval;
    }
    
    
    
    
    virtual mpfr FreshEval(mpfr) override
    {
        mpfr retval{1};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_mult_or_div_[ii])
            {
                retval *= children_[ii]->Eval<mpfr>();
            }
            else
            {
                retval /= children_[ii]->Eval<mpfr>();
            }
        }
        
        return retval;
    }
    
    
    
    
    
    
private:
    // Stores the mult/div of a factor.  There is a one-one
    // correspondence between elements of children_sign_ and children_.  This
    // is enforced by the AddChild method, redefined in MultOperator.
    
    // TODO(JBC): If we add method to delete child, must also delete children_mult_ entry.
    std::vector<bool> children_mult_or_div_;

};


} // re: namespace bertini



using Node = bertini::Node;
using MultOperator = bertini::MultOperator;

inline std::shared_ptr<Node> operator*=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
{
	lhs->AddChild(rhs);
	return lhs;
}


inline std::shared_ptr<Node> operator/=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
{
	lhs->AddChild(rhs,false);
	return lhs;
}

inline std::shared_ptr<Node>& operator*=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
{
	
	std::shared_ptr<Node> temp = std::make_shared<MultOperator>();
	
	std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
	std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(rhs);
	
	lhs.swap(temp);
	return lhs;
}


inline std::shared_ptr<Node>& operator/=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
{
	
	std::shared_ptr<Node> temp = std::make_shared<MultOperator>();
	
	std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
	std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(rhs,false);
	
	lhs.swap(temp);
	return lhs;
}


inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
{
	return lhs*=rhs;
}


inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
{
	return lhs/=rhs;
}




#endif
