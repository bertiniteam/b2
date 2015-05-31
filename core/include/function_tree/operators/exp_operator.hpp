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
// exp_operator.hpp:  Declares the class  ExpOperator.



#ifndef b2Test_ExpOperator_h
#define b2Test_ExpOperator_h

#include "function_tree/Node.hpp"
#include "function_tree/operators/unary_operator.hpp"



namespace bertini {

	using std::pow;
	
	/**
	 Node -> UnaryOperator -> ExpOperator

	 Description: This class represents the exponentiation operator.  The base is stored in
	 children_, and an extra variable(exponent_) stores the exponent.  FreshEval is
	 defined as the exponention operation.
	 */
class ExpOperator : public virtual UnaryOperator
{
public:
    // These do nothing for a constant
    std::string PrintNode() override
    {
        return child_->PrintNode() + "^" + std::to_string(exponent());
    }
	
	
	
	/**
	 Virtual polymorphic method for printing to an arbitrary stream.
	 */
	virtual void print(std::ostream & target) const override
	{
		target << "(";
		child_->print(target);
		target << "^" << exponent() << ")";
	}
		
	/**
	 Get the integet exponent of an ExpOperator
	 */
    void set_exponent(int exp)
    {
        exponent_ = exp;
    }
	
	
	/**
	 Get the exponent of an ExpOperator
	 */
    int exponent() const
    {
        return exponent_;
    }

	virtual ~ExpOperator() = default;

	
	/**
	 Constructor, passing in the Node you want as the base, and the integer you want for the power.
	 */
	ExpOperator(const std::shared_ptr<Node> & N, int p = 1) : exponent_(p), UnaryOperator(N)
	{}
	
	
	ExpOperator(){}

protected:
    // Specific implementation of FreshEval for exponentiate.
    // TODO(JBC): How do we implement exp for more complicated types?
    dbl FreshEval(dbl) override
    {
        return pow(child_->Eval<dbl>(), exponent_);
    }
    
    mpfr FreshEval(mpfr) override
    {
        return pow(child_->Eval<mpfr>(),exponent_);
    }

    
    
    
    
    
private:
	
	int exponent_ = 1; ///< Exponent for the exponenetial operator
};


}// re: namespace bertini



namespace {

	using Node = bertini::Node;
	using ExpOperator = bertini::ExpOperator;

	/**
	 Overloading of the power function to create an ExpOperator with base N and power p.
	 */
	inline std::shared_ptr<Node> pow(const std::shared_ptr<Node> & base, int power)
	{
		return std::make_shared<ExpOperator>(base,power);
	}

}
#endif
