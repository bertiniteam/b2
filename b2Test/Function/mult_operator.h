// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// mult_operator.h:  Declares the class MultOperator.



#ifndef b2Test_MultOperator_h
#define b2Test_MultOperator_h

#include "node.h"
#include "nary_operator.h"



// Node -> NaryOperator -> MultOperator
// Description: This class represents multiplication operator.  All children are factors and are stored
// in a vector.  FreshEval method is defined for multiplication.
class MultOperator : public NaryOperator
{
public:
    // See node.h for description
    // TODO(JBC): Implement this
    std::string PrintNode() override {return "";}






protected:
    // Specific implementation of FreshEval for multiply.
    dbl FreshEval(dbl) override
    {
        dbl retval{1};
        for(auto& vv : children_)
        {
            retval *= vv->Eval<dbl>();
        }
        
        return retval;
    }
    
    mpfr FreshEval(mpfr) override
    {
        mpfr retval{1};
        for(auto& vv : children_)
        {
            retval *= vv->Eval<mpfr>();
        }
        
        return retval;
    }
};




#endif
