//
//  Operand.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__Operand__
#define __Bertini2__Operand__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class Operand
{
protected:
    bool isLeaf;
    bool isExp;

    
public:
    int IsLeaf() { return isLeaf;};
    int IsExp() {return isExp;};
    virtual double evaluate() = 0;
    virtual double evaluate(double* vars) = 0;
    virtual std::stringstream print() = 0;
    virtual std::vector<Operand*> getOperands() = 0;
    virtual int getNumOperands() = 0;
    virtual void addOperand(Operand* op) = 0;
};


#endif /* defined(__Bertini2__Operand__) */
