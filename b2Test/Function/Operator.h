//
//  Operator.h
//  b2Test
//
//  Created by Collins, James B. on 5/2/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_Operator_h
#define b2Test_Operator_h

#include "Node.h"

class Operator : public Node
{
protected:
    std::tuple< std::pair<dbl,bool>, std::pair<mpfr,bool> > current_value;
    
public:
    // Delete add_child methods
    virtual void add_Child(std::unique_ptr<Node>) = 0;
    

};

#endif
