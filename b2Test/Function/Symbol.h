//
//  Symbol.h
//  b2Test
//
//  Created by Collins, James B. on 5/1/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef b2Test_Symbol_h
#define b2Test_Symbol_h


#include "Node.h"

class Symbol : public Node
{
   //I forgot why we needed this class
public:
    virtual void add_Child(std::unique_ptr<Node> child) {};
    virtual void add_Child(Node* child) {};


    

};



#endif
