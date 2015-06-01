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
// constant.hpp:  Declares the class Number.



#ifndef b2Test_Number_h
#define b2Test_Number_h


#include "function_tree/symbols/symbol.hpp"


namespace bertini {



// Node -> Symbol -> Number
// Description: This class represents constant leaves to a function tree.  FreshEval simply returns
// the value of the constant.
class Number : public virtual Node
{
public:
    Number(){};
	
	
    Number(double val)
    {
        std::get< std::pair<dbl,bool> >(current_value_).first = val;
        std::get< std::pair<mpfr,bool> >(current_value_).first = val;
    }
	
	
    // Ctor that reads in a string and converts it to a number of the appropriate type
    Number(std::string sval)
    {
        std::get< std::pair<dbl,bool> >(current_value_).first = stod(sval);
        std::get< std::pair<mpfr,bool> >(current_value_).first = stod(sval);
    }
    
	
	
	
	
	
	
    // These do nothing for a constant
    std::string PrintNode() override {return "";}

	virtual void print(std::ostream & target) const override
	{
		target << std::get< std::pair<mpfr,bool> >(current_value_).first;
	}
	
	virtual void Reset() override
	{
		// nothing to reset here
	}
	
	
	virtual ~Number() = default;

	
protected:
    // Return value of constant
    dbl FreshEval(dbl) override
    {
        return std::get< std::pair<dbl,bool> >(current_value_).first;
    }
    
    mpfr FreshEval(mpfr) override
    {
        return std::get< std::pair<mpfr,bool> >(current_value_).first;
    }
};
	
	
} // re: namespace bertini

#endif
