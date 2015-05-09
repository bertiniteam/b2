// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// exp_operator.h:  Declares the class  ExpOperator.



#ifndef b2Test_ExpOperator_h
#define b2Test_ExpOperator_h

#include "Node.h"
#include "unary_operator.h"



// Node -> UnaryOperator -> ExpOperator
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
