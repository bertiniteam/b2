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
// variable.h:  Declares the class Variable.



#ifndef b2Test_Variable_h
#define b2Test_Variable_h

#include "symbol.h"




// Node -> Symbol -> Variable
// Description: This class represents variable leaves in the function tree.  FreshEval returns
// the current value of the variable.
class Variable : public Symbol
{
public:
    Variable(){};
    Variable(double val)
    {
        std::get< std::pair<dbl,bool> >(current_value_).first = val;
        std::get< std::pair<mpfr,bool> >(current_value_).first = val;
    }
    
    
    // These do nothing for a constant
    std::string PrintNode() override {return "";}
    
    
    
    
    
    // This sets the value for the variable
    template <typename T>
    void set_current_value(T val)
    {
        std::get< std::pair<T,bool> >(current_value_).first = val;
        std::get< std::pair<T,bool> >(current_value_).second = false;
    }







protected:
    // Return current value of the variable.
    dbl FreshEval(dbl) override
    {
        return std::get< std::pair<dbl,bool> >(current_value_).first;
    }
    
    mpfr FreshEval(mpfr) override
    {
        return std::get< std::pair<mpfr,bool> >(current_value_).first;
    }
};


#endif
