//
//  SymbOp.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__SymbOp__
#define __Bertini2__SymbOp__

#include <iostream>
#include <vector>
#include "Operand.h"


class OpOperand : public Operand
{
protected:
    std::vector<Operand*> operands;
    int numOperands;
    
public:
    virtual double evaluate() = 0;
    virtual double evaluate(double* vars) = 0;
    virtual std::stringstream print() = 0;
    
    int getNumOperands() {return numOperands;};
    virtual std::vector<Operand*> getOperands() {return operands;};
    
    virtual void addOperand(Operand* op) {if(op != 0) {operands.push_back(op); numOperands++;}};
};




#endif /* defined(__Bertini2__SymbOp__) */
