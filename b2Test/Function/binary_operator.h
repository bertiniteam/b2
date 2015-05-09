// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// binary_operator.h:  Declares the class BinaryOperator.



#ifndef b2Test_BinaryOperator_h
#define b2Test_BinaryOperator_h

#include <vector>

#include "node.h"



// Description: This class is an interface for all binary operators, such as division.
// Children of the operator are stored in a vector.
class BinaryOperator : public Node
{
public:
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        if(children.size() == 2)
        {
            std::cout << "Already have two children on a binary node, can't add more\n But I don't know how to handle errors yet(JBC)\n";
            return;
        }
        children.push_back(std::move(child));
    }




protected:
    //Stores the two children for this binary node
    //A vector is used so that it is easy to add a
    //child.
    //TODO(JBC): Should this be changed for speed sake?
    std::vector< std::unique_ptr<Node> > children_;
    
    
    
    
    
    
    
    ////////////// TESTING /////////////////
    virtual void PrintTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
        tabcount++;
        for(auto vv : children)
        {
            vv->PrintTree();
        }
        tabcount--;
    }
    ////////////// TESTING /////////////////
    
};




#endif
