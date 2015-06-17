//This file is part of Bertini 2.0.
//
//function.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 6/11/15.
//
//
// jacobian.hpp:  Declares the class Jacobian.


#ifndef b2Test_Jacobian_h
#define b2Test_Jacobian_h


#include "function_tree/node.hpp"



namespace bertini {
    
    
    /**
     Node -> Function
     This class defines a function.  It stores the entry node for a particular functions tree.
     */
    class Jacobian : public Function
    {
    public:
        
        
        /**
         The default constructor
         */
        Jacobian()
        {};
        
        
        /**
         Constructor defines entry node at construct time.
         */
        Jacobian(const std::shared_ptr<Node> & entry) : Function(entry)
        {
        }
        
        
       
        
        
        
        
        virtual ~Jacobian() = default;
        
        
        
        
    };
    
    
} // re: namespace bertini



#endif
