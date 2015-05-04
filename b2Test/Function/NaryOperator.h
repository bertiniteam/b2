//
//  NaryOperator.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_NaryOperator_h
#define b2Test_NaryOperator_h

#include <vector>

#include "Node.h"


class NaryOperator : public Node
{
protected:
    std::vector< std::unique_ptr<Node> > children;
    
    //Deleted functions
    
    
public:
    virtual void add_Child(std::unique_ptr<Node> child) override
    {
        children.push_back(std::move(child));
    }
    
    virtual void add_Child(Node* child) override
    {
        children.push_back( std::move( std::unique_ptr<Node>(child) ) );
    }

};




#endif
