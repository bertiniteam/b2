//
//  ExpOp.cpp
//  Bertini2
//
//  Created by James Collins on 11/20/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include <cmath>
#include "ExpOp.h"
#include "LeafOperand.h"


ExpOp::ExpOp(Operand* inBase, int inExp)
{
    base = inBase;
    exp = inExp;
    type = TYPE_EXP;
}




//double ExpOp::evaluate(double* vars)
//{
//    double retval = pow(base->evaluate(vars),exp);
//    
//    return retval;
//}


Symbol* ExpOp::evaluate(Symbol* vars[])
{
    Symbol* retval = base->evaluate(vars)->clone();
    
    retval = retval->exp(exp);
    
    return retval;
  
}


std::stringstream ExpOp::print()
{
    std::stringstream ret;
    std::stringstream term = base->print();
    
    if(base->getType() == TYPE_LEAF)
    {
        ret  << term.str()  << "^" << exp;
    }
    else
    {
        ret  << "(" << term.str()  << ")" << "^" << exp;
    }
    
        
    return ret;
}





Operand* ExpOp::diff(int varIndex)
{
    if(exp == 0)
    {
        return new LeafOperand(new IntSymb(0));
    }
    else if(exp == 1)
    {
        return base->diff(varIndex);
    }
    else
    {
        MultOp* retMult = new MultOp();
        retMult->addOperand(new LeafOperand(new IntSymb(exp)));
        retMult->addOperand(base->diff(varIndex));
        retMult->addOperand(new ExpOp(base,exp-1));
        return retMult;
    }
    
    
    
}







/***************  cleanMults()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Addition node.  Recursively calls each branch to clean
 mults for nodes lower down the tree.
 ****************************************************/

void ExpOp::cleanMults()
{
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->cleanMults();
        if(operands[ii]->getType() == TYPE_MULT)
        {
            if(operands[ii]->getNumOperands() == 1)
            {
                if(operands[ii]->getOperands()[0]->getType() == TYPE_LEAF)
                {
                    LeafOperand* leaf = (LeafOperand*)operands[ii]->getOperands()[0];
                    if(leaf->getIsZero())
                    {
                        delete operands[ii];
                        operands[ii] = new LeafOperand(0,1);
                    }
                }
            }
        }
    }
}


