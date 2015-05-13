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
