// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// symbol.h:  Declares the class Symbol.

#ifndef b2Test_Symbol_h
#define b2Test_Symbol_h


#include "node.h"


// Description: This class could work as an interface for all non-operators.
// NOTE: Currently this does nothing!!!!
class Symbol : public Node
{
   //I forgot why we needed this class
public:
    virtual void AddChild(std::shared_ptr<Node> child) override {};


    

};



#endif
