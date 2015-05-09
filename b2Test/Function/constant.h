// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// constant.h:  Declares the class  Constant.



#ifndef b2Test_Constant_h
#define b2Test_Constant_h


#include "symbol.h"



// Node -> Symbol -> Constant
// Description: This class represents constant leaves to a function tree.  FreshEval simply returns
// the value of the constant.
class Constant : public Symbol
{
public:
    Constant(){};
    Constant(double val)
    {
        std::get< std::pair<dbl,bool> >(current_value_).first = val;
        std::get< std::pair<mpfr,bool> >(current_value_).first = val;
    }
    
    
    // These do nothing for a constant
    std::string PrintNode() override {return "";}








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

#endif
