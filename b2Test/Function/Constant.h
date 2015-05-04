//
//  Constant.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_Constant_h
#define b2Test_Constant_h


#include "Node.h"
#include "Symbol.h"


class Constant : public Symbol
{
protected:
    dbl fresh_eval(dbl) override
    {
        return std::get< std::pair<dbl,bool> >(current_value).first;
    }
    
    mpfr fresh_eval(mpfr) override
    {
        return std::get< std::pair<mpfr,bool> >(current_value).first;
    }
    
    
    

    
public:
    
    // These do nothing for a constant
    std::string get_string() override {return "";}
    
};

#endif
