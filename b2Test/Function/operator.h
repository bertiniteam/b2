// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// operator.h:  Declares the class Operator.

#ifndef b2Test_Operator_h
#define b2Test_Operator_h

#include "node.h"


// Description: This class is could work as an interface for all operators in a function tree.
// NOTE: Currently this class is not in use!!!!!
class Operator : public Node
{
    
public:
    // Delete add_child methods
    virtual void AddChild(std::shared_ptr<Node>) = 0;
    

};

#endif
