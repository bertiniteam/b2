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
// exp_operator.h:  Declares the class  ExpOperator.



#ifndef b2Test_ExpOperator_h
#define b2Test_ExpOperator_h

#include "Node.h"
#include "unary_operator.h"



// Node -> UnaryOperator -> ExpOperator
//
// Description: This class represents the exponentiation operator.  The base is stored in
// children_, and an extra variable(exponent_) stores the exponent.  FreshEval is
// defined as the exponention operation.
class ExpOperator : public UnaryOperator
{
public:
    // These do nothing for a constant
    std::string PrintNode() override
    {
        return children_->PrintNode() + "^" + std::to_string(exponent());
    }
    
    
    void set_exponent(int exp)
    {
        exponent_ = exp;
    }
    
    int exponent()
    {
        return exponent_;
    }





protected:
    // Specific implementation of FreshEval for exponentiate.
    // TODO(JBC): How do we implement exp for more complicated types?
    dbl FreshEval(dbl) override
    {
        return pow(children_->Eval<dbl>(), exponent_);
    }
    
    mpfr FreshEval(mpfr) override
    {
        return pow(children_->Eval<mpfr>(),exponent_);
    }

    
    
    
    
    
private:
    // Exponent for the exponenetial operator
    int exponent_ = 1;
};




#endif
