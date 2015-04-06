//
//  AdditionOp.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__AdditionOp__
#define __Bertini2__AdditionOp__

#include <iostream>
#include <vector>
#include "OpOperand.h"


class AdditionOp : public OpOperand
{
public:
    AdditionOp() {operands = std::vector<Operand*>();};
    AdditionOp(int nOperands, std::vector<Operand*> inOperands);
    virtual double evaluate() {return 1;};
    virtual double evaluate(double* vars);
    virtual std::stringstream print();
    
    
};



#endif /* defined(__Bertini2__AdditionOp__) */
