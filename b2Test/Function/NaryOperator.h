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
    std::vector< std::shared_ptr<Node> > children;
    
    
    
public:
    // Add a child onto the container for this operator
    virtual void add_Child(std::shared_ptr<Node> child) override
    {
        children.push_back(std::move(child));
    }
    
    
    int num_Children()
    {
        return children.size();
    }
    
    std::shared_ptr<Node> get_first_child()
    {
        return children[0];
    }
    
    
    ////////////// TESTING /////////////////
    virtual void printTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->eval<dbl>() << std::endl;
        tabcount++;
        for(auto vv : children)
        {
            vv->printTree();
        }
        tabcount--;
    }
    ////////////// TESTING /////////////////

};




#endif
