//
//  DivideOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_DivideOperator_h
#define b2Test_DivideOperator_h

#include "Node.h"
#include "binary_operator.h"


class DivideOperator : public BinaryOperator
{
public:
    // These do nothing for a constant
    std::string PrintNode() override
    {
        return "(" + children_[0].PrintNode() + "/" + children_[1].PrintNode();
    }





protected:
    // Specific implementation of FreshEval for divide.
    dbl FreshEval(dbl) override
    {
        return children_[0].Eval<dbl>()/children_[1].Eval<dbl>();
    }
    
    mpfr FreshEval(mpfr) override
    {
        return children_[0].Eval<mpfr>()/children_[1].Eval<mpfr>();
    }
};




#endif
