//
//  NegateOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_NegateOperator_h
#define b2Test_NegateOperator_h

#include "Node.h"
#include "UnaryOperator.h"


class NegateOperator : public UnaryOperator
{
protected:
    dbl fresh_eval(dbl) override
    {
        return (-1)*children->eval<dbl>();
    }
    
    mpfr fresh_eval(mpfr) override
    {
        return (-1)*children->eval<mpfr>();
    }
    
    
    
    
    
public:
    // These do nothing for a constant
    std::string get_string() override
    {
        return "-" + children->get_string();
    }
};



#endif
