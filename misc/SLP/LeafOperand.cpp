//
//  LeafOperand.cpp
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "LeafOperand.h"


double LeafOperand::evaluate(double* vars)
{
    if(isVar)
    {
        return vars[var];
    }
    else
    {
        return value;
    }
    
};


std::stringstream LeafOperand::print()
{
    std::stringstream ret;
    
    if(isVar)
    {
        ret << "x" << var;
    }
    else
    {
        ret << value;
    }
    
    return ret;
};