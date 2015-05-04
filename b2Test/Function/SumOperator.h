//
//  SumOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_SumOperator_h
#define b2Test_SumOperator_h

#include <vector>
#include <string>   

#include "Node.h"
#include "NaryOperator.h"


class SumOperator : public NaryOperator
{
private:
    std::vector<bool> children_sign;
    
protected:
    std::tuple< std::pair<dbl,bool>, std::pair<mpfr,bool> > current_value;
    
    
    dbl fresh_eval(dbl) override
    {
        dbl retval{};
        for(auto& vv : children)
        {
            retval += vv->eval<dbl>();
        }
        
        return retval;
    }
    
    mpfr fresh_eval(mpfr) override
    {
        mpfr retval{};
        for(auto& vv : children)
        {
            retval += vv->eval<mpfr>();
        }
        
        return retval;
    }
    
    void add_Child(std::unique_ptr<Node> child, bool sign)
    {
        children.push_back(std::move(child));
    }
    
    
    
    
public:
    // These do nothing for a constant
    std::string get_string() override {return "";}
};


#endif
