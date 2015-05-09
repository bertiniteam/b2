// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// negate_operator.h:  Declares the class  NegateOperator.



#ifndef b2Test_NegateOperator_h
#define b2Test_NegateOperator_h

#include "Node.h"
#include "unary_operator.h"



// Node -> UnaryOperator -> NegateOperator
// Description: This class represents the negation operator.  FreshEval method
// is defined for negation and multiplies the value by -1.
class NegateOperator : public UnaryOperator
{
public:
    // These do nothing for a constant
    std::string PrintNode() override
    {
        return "-" + children_->PrintNode();
    }






protected:
    // Specific implementation of FreshEval for negate.
    dbl FreshEval(dbl) override
    {
        return (-1)*children_->Eval<dbl>();
    }
    
    mpfr FreshEval(mpfr) override
    {
        return (-1)*children_->Eval<mpfr>();
    }
};
#endif
