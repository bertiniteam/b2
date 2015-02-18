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
    isExp = true;
}




double ExpOp::evaluate(double* vars)
{
    double retval = pow(base->evaluate(vars),exp);
    
    return retval;
}


std::stringstream ExpOp::print()
{
    std::stringstream ret;
    std::stringstream term = base->print();
    
    ret << "(" << term.str() << ")" << "^" << exp;
        
    return ret;
}