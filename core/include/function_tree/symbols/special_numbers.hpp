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
// special_numbers.hpp:  Declares the class SpecialNumber.



#ifndef b2Test_Number_h
#define b2Test_Number_h

#include <boost/math/constants/constants.hpp>

#include "function_tree/symbols/symbol.hpp"


namespace bertini {
    
    
    
    // Node -> Symbol -> Number
    // Description: This class represents constant leaves to a function tree.  FreshEval simply returns
    // the value of the constant.
    class SpecialNumber : public virtual NamedSymbol
    {
    public:
        SpecialNumber(){};
        
        
        /**
         Input string num denotes the type of special number this node is.
         "pi" or "Pi" is the number pi
         "e" or "E" is exp(1)
         "i" or "I" or "1i" is the imaginary number
         */
        SpecialNumber(std::string num)
        {
            using namespace boost::math::constants;
            
            switch (num)
            {
                case "pi":
                case "Pi":
                    std::get< std::pair<dbl,bool> >(current_value_).first = dbl(pi<dbl>,0.0);
                    std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(pi<mpfr>,0.0);
                    name("pi");
                    break;
                case "e":
                case "E":
                    std::get< std::pair<dbl,bool> >(current_value_).first = dbl(e<dbl>,0.0);
                    std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(e<mpfr>,0.0);
                    name("e");
                    break;
                case "i":
                case "I":
                case "1i":
                    std::get< std::pair<dbl,bool> >(current_value_).first = dbl(0.0,1.0);
                    std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(0.0,1.0);
                    name("i");
                    break;
            }
        }
        
        
        
        
        
        
        
    
    virtual void print(std::ostream & target) const override
    {
    target << std::get< std::pair<mpfr,bool> >(current_value_).first;
    }
    
    virtual void Reset() override
    {
        // nothing to reset here
    }
    
    
    virtual ~SpecialNumber() = default;
    
    
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
