//
//  BinaryOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_BinaryOperator_h
#define b2Test_BinaryOperator_h

#include <vector>

#include "Operator.h"


class BinaryOperator : public Operator
{
protected:
    std::vector< std::unique_ptr<Node> > children;
    
public:
    virtual void add_Child(std::shared_ptr<Node> child) override
    {
        if(children.size() == 2)
        {
            std::cout << "Already have two children on a binary node, can't add more\n But I don't know how to handle errors yet(JBC)\n";
            return;
        }
        children.push_back(std::move(child));
    }
    
};




#endif
