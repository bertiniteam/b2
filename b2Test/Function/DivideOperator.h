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
#include "BinaryOperator.h"


class DivideOperator : public BinaryOperator
{
protected:
    dbl fresh_eval(dbl) override
    {
        return children[0].eval<dbl>()/children[1].eval<dbl>();
    }
    
    mpfr fresh_eval(mpfr) override
    {
        return children[0].eval<mpfr>()/children[1].eval<mpfr>();
    }
    
    
    
    
    
public:
    // These do nothing for a constant
    std::string get_string() override
    {
        return "(" + children[0].get_string() + "/" + children[1].get_string();
    }
};




#endif
