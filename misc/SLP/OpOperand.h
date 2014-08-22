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
#include "Symbol.h"


class OpOperand : public Operand
{
protected:
    std::vector<Operand*> operands;
    int numOperands;
    
public:
    virtual Symbol* evaluate(Symbol* vars[]) = 0;
    virtual std::stringstream print() = 0;
    virtual Operand* diff(int varIndex) = 0;
    
    int getNumOperands() {return numOperands;};
    virtual std::vector<Operand*> getOperands() {return operands;};
    
    virtual void addOperand(Operand* op) {if(op != 0) {operands.push_back(op); numOperands++;}};
    
    
    void clear() {operands.clear(); numOperands = 0;};
    
    virtual void combineAdds() = 0;
    virtual void combineMults() = 0;
    virtual void cleanMults() = 0;
};




#endif /* defined(__Bertini2__SymbOp__) */
