// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// nary_operator.h:  Declares the class NaryOperator.

#ifndef b2Test_NaryOperator_h
#define b2Test_NaryOperator_h

#include <vector>

#include "node.h"




// Description: This class is an interface for all n-ary operators, such as summation and multiplication.
// Children of the operator are stored in a vector and methods to add and access children are available
// in this interface.
class NaryOperator : public Node
{
public:
    // Add a child onto the container for this operator
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        children_.push_back(std::move(child));
    }
    
    
    int children_size()
    {
        return children_.size();
    }
    
    std::shared_ptr<Node> first_child()
    {
        return children_[0];
    }
    
    
    ////////////// TESTING /////////////////
    virtual void PrintTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << tabcount+1 << "." <<  boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>() << std::endl;
        tabcount++;
        for(auto vv : children_)
        {
            vv->PrintTree();
        }
        tabcount--;
    }
    ////////////// TESTING /////////////////
    
    
    
    
    
protected:
    //Stores all children for this operator node.
    //This is an NaryOperator and can have any number of children.
    std::vector< std::shared_ptr<Node> > children_;
    
    
    

};




#endif
