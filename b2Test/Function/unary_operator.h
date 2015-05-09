// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// unary_operator.h:  Declares the class UnaryOperator.



#ifndef b2Test_UnaryOperator_h
#define b2Test_UnaryOperator_h


#include <vector>

#include "node.h"


// Description: This class is an interface for all unary operators, such as negation.
// The sole child is stored in a shared_ptr.
class UnaryOperator : public Node
{
public:
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        children_ = std::move(child);
    }
    
    
    //Return the only child for the unary operator
    std::shared_ptr<Node> first_child()
    {
        return children_;
    }
    
    
    ////////////// TESTING /////////////////
    virtual void PrintTree()
    {
        for(int ii = 0; ii < tabcount; ++ii)
        {
            std::cout << "\t";
        }
        std::cout << boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
        tabcount++;
        children_->PrintTree();
        tabcount--;
    }
    ////////////// TESTING /////////////////






protected:
    //Stores the single child of the unary operator
    std::shared_ptr<Node> children_;
};



#endif
