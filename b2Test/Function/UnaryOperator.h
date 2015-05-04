//
//  UnaryOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_UnaryOperator_h
#define b2Test_UnaryOperator_h


#include <vector>

#include "Operator.h"


class UnaryOperator : public Operator
{
protected:
    std::unique_ptr<Node> children;
    
public:
    virtual void add_Child(std::unique_ptr<Node> child) override
    {
        children = std::move(child);
    }
    

};



#endif
