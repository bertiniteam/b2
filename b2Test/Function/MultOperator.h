//
//  MultOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_MultOperator_h
#define b2Test_MultOperator_h

#include "Node.h"
#include "NaryOperator.h"


class MultOperator : public NaryOperator
    {
        
    protected:
        
        
        dbl fresh_eval(dbl) override
        {
            dbl retval{1};
            for(auto& vv : children)
            {
                retval *= vv->eval<dbl>();
            }
            
            return retval;
        }
        
        mpfr fresh_eval(mpfr) override
        {
            mpfr retval{1};
            for(auto& vv : children)
            {
                retval *= vv->eval<mpfr>();
            }
            
            return retval;
        }
        
        
        
        
        
    public:
        // These do nothing for a constant
        std::string get_string() override {return "";}
    };




#endif
