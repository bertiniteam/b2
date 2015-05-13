//This file is part of Bertini 2.0.
//
//Foobar is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Foobar is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
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
        std::cout << tabcount+1 << "." <<  boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
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
