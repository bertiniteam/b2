//
//  LeafOperand.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__LeafOperand__
#define __Bertini2__LeafOperand__

#include <iostream>
#include "Operand.h"


class LeafOperand : public Operand
{
private:
    double value;
    int var;
    int isVar;
    
public:
    LeafOperand() {isLeaf = 1; isVar = 0; isExp = 0;};
    LeafOperand(double inval) {value = inval; isLeaf = 1; isVar = 0; var = -1;};
    LeafOperand(int inVar) {var = inVar; isLeaf = 1; isVar = 1; value = 0.0;};
    
    virtual double evaluate() { return value;};
    virtual double evaluate(double* vars);
    void setValue(double inval) {value = inval; isVar = 0;};
    void setVar(int inVar) {var = inVar; isVar = 1;};
    
    virtual std::stringstream print();
    
    virtual std::vector<Operand*> getOperands() {std::vector<Operand*> ret; return ret;};
    virtual int getNumOperands() {return 0;}
    
    virtual void addOperand(Operand* op) {};
};

#endif /* defined(__Bertini2__LeafOperand__) */
