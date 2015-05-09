// BoilerPlate licence
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
