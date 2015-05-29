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
// negate_operator.h:  Declares the class  NegateOperator.



#ifndef b2Test_NegateOperator_h
#define b2Test_NegateOperator_h

#include "function_tree/Node.h"
#include "function_tree/operators/unary_operator.h"


namespace bertini {

// Node -> UnaryOperator -> NegateOperator
// Description: This class represents the negation operator.  FreshEval method
// is defined for negation and multiplies the value by -1.
class NegateOperator : public  virtual UnaryOperator
{
public:
	
	NegateOperator(){}
	
	NegateOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
	{};
	
    // These do nothing for a constant
    std::string PrintNode() override
    {
        return "-" + child_->PrintNode();
    }



	virtual void print(std::ostream & target) const override
	{
		target << "-(";
		child_->print(target);
		target << ")";
	}

	virtual ~NegateOperator() = default;

protected:
    // Specific implementation of FreshEval for negate.
    dbl FreshEval(dbl) override
    {
        return (-1)*child_->Eval<dbl>();
    }
    
    mpfr FreshEval(mpfr) override
    {
        return (-1)*child_->Eval<mpfr>();
    }
};

} // re: namespace bertini


namespace  {

	using Node = bertini::Node;
	using NegateOperator = bertini::NegateOperator;
	
inline std::shared_ptr<Node> operator-(const std::shared_ptr<Node> & rhs)
{
	return std::make_shared<NegateOperator>(rhs);
}

}

#endif
