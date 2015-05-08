//
//  ExpOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/7/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_ExpOperator_h
#define b2Test_ExpOperator_h

#include "Node.h"
#include "UnaryOperator.h"


class ExpOperator : public UnaryOperator
{
private:
    int exponent = 1;
    
protected:
    dbl fresh_eval(dbl) override
    {
        return pow(children->eval<dbl>(), exponent);
    }
    
    mpfr fresh_eval(mpfr) override
    {
        return pow(children->eval<mpfr>(),exponent);
    }
    
    
    
    
    
public:
    // These do nothing for a constant
    std::string get_string() override
    {
        return "-" + children->get_string();
    }
    
    
    void setExponent(int exp)
    {
        exponent = exp;
    }
};




#endif
