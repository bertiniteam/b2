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
// function.h:  Declares the class Function.


#ifndef b2Test_Function_h
#define b2Test_Function_h


#include "node.h"




// Node -> Function
// This class defines a function.  It stores the entry node for a particular functions tree.
// TODO(JBC): What else will be in this class?
class Function : public Node
{
public:
    // Constructor defines entry node
    Function(std::shared_ptr<Node> entry)
    {
        entry_node_ = std::move(entry);
    }
    
    
    // This returns a string that holds a print out of the function.
    virtual std::string PrintNode()
    {
        return entry_node_->PrintNode();
    }
    
    // Add a child onto the container for this operator
    virtual void AddChild(std::shared_ptr<Node> entry) override
    {
        entry_node_ = entry;
    }
    
    
    // Add a child onto the container for this operator
    void AddChild(Node* entry)
    {
        entry_node_ = std::shared_ptr<Node>(entry);
    }

    
    
    std::shared_ptr<Node> entry_node() {return entry_node_;};

    
    
    ////////////// TESTING /////////////////
    virtual void PrintTree()
    {
        entry_node_->PrintTree();
    }
    ////////////// TESTING /////////////////

private:
    // Calls FreshEval on the entry node to the tree.
    virtual dbl FreshEval(dbl d)
    {
        return entry_node_->Eval<dbl>();
    }
    
    virtual mpfr FreshEval(mpfr m)
    {
        return entry_node_->Eval<mpfr>();
    }
    
    
    // The top node for the function.  This should be a SumOperator
    std::shared_ptr<Node> entry_node_;
};






#endif
