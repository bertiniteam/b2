//
//  LeafOperand.cpp
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "LeafOperand.h"



/***************  evaluate(Symbol* vars[])  ********************
 Input: Array to Symbol pointers representing values of variables
 
 Output: Symbol pointer to the value of the leaf
 
 Desc: Evaluates the leaf node.
 **********************************************/

Symbol* LeafOperand::evaluate(Symbol* vars[])
{
    Symbol* retval;
    if(isVar)
    {
        retval = vars[varIndex]->clone();
    }
    else
    {
        retval = symb->clone();
    }
    
    return retval;
};











/***************  print()  ********************
 Input: None
 
 Output: stringstream containing character output of polynomial
 
 Desc: Prints the Leaf node.  If a variable, prints the variable name,
 if a coefficient, prints the number.
 **********************************************/

std::stringstream LeafOperand::print()
{
    std::stringstream ret;
    
    
    if(isVar)
    {
        ret << "x" << varIndex;
    }
    else
    {
        ret << symb->print().str() << "*";
    }
    
    return ret;
};









/***************  diff(int varIndex)  ********************
 Input: the index of the variable to be differentiated with respect to
 
 Output: Leaf node with differentiated Leaf
 
 Desc: Differentiates the Leaf node.  Is either a 1 or 0, depending on
 whether leaf is a variable or coefficient.
 **********************************************/

LeafOperand* LeafOperand::diff(int invarIndex)
{
    if(isVar)
    {
        if(varIndex == invarIndex)
        {
            return new LeafOperand(1,0);
        }
        else
        {
            return new LeafOperand(0,1);
        }
    }
    else
    {
        return new LeafOperand(0,1);
    }
     
}







