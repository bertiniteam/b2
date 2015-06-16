//This file is part of Bertini 2.0.
//
//number.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//number.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with number.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// number.hpp:  Declares the class Number.



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
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(val,0.0);
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(val,0.0);
		}
        
        Number(double rval, double ival)
        {
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(rval,ival);
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(rval,ival);
		}
		
		
		
		// Ctor that reads in a string for the real component and converts it to a number of the appropriate type
		Number(std::string sval)
		{
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(stod(sval), 0.0);
			std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(sval, "0.0");
		}
		
        // Ctor that reads in two strings for a complex number and converts it to a number of the appropriate type
        Number(std::string srval, std::string sival)
        {
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(stod(srval), stod(sival));
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(srval, sival);
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
        
        
        /**
         Differentiates a number.  Should this return the special number Zero?
         */
        virtual std::shared_ptr<Node> Differentiate() const override
        {
            return std::make_shared<Number>(0.0);
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
