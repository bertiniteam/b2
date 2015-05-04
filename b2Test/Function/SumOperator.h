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
    
    
    dbl fresh_eval(dbl) override
    {
        dbl retval{};
        for(int ii = 0; ii < children.size(); ++ii)
        {
            if(children_sign[ii])
            {
                retval += children[ii]->eval<dbl>();
            }
            else
            {
                retval -= children[ii]->eval<dbl>();
            }
        }
        
        return retval;
    }
    
    mpfr fresh_eval(mpfr) override
    {
        mpfr retval{};
        for(int ii = 0; ii < children.size(); ++ii)
        {
            if(children_sign[ii])
            {
                retval += children[ii]->eval<mpfr>();
            }
            else
            {
                retval -= children[ii]->eval<mpfr>();
            }
        }
        
        return retval;
    }
    
    
    
    
    
public:
    // These do nothing for a constant
    std::string get_string() override {return "";}
    
    virtual void add_Child(std::unique_ptr<Node> child) override
    {
        children.push_back(std::move(child));
        children_sign.push_back(true);
    }

    void add_Child(std::unique_ptr<Node> child, bool sign)
    {
        children.push_back(std::move(child));
        children_sign.push_back(sign);
    }
};


#endif
