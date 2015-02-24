//
//  MultOp.cpp
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "MultOp.h"
#include "LeafOperand.h"


MultOp::MultOp(int nOperands, std::vector<Operand*> inOperands)
{
    numOperands = nOperands;
    operands = inOperands;
    isExp = false;
}




double MultOp::evaluate(double* vars)
{
    double retval = 1;

    for(int ii = 0; ii < numOperands; ii++)
    {

        if(operands[ii]->IsLeaf())
        {
            retval *= operands[ii]->evaluate(vars);
        }
        else
        {
            retval *= operands[ii]->evaluate(vars);
        }
    }
    
    return retval;
}



std::stringstream MultOp::print()
{
    std::stringstream ret;
    if(operands.size() > 0)
    {
        std::stringstream term = operands[0]->print();
        
        ret << "(" << term.str() << ")";
        
        for(int ii = 1; ii < operands.size(); ii++)
        {
            term = operands[ii]->print();
            ret << "(" << term.str() << ")";
        }
    }
    
    return ret;
}