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
#include "Operand.h"


class AdditionOp : public Operand
{
private:
    std::vector<int> coeffs;
    
public:
    AdditionOp() {operands = std::vector<Operand*>(); numOperands = 0; type = TYPE_ADD;};
    AdditionOp(std::vector<Operand*> inOperands);
    AdditionOp(std::vector<Operand*> inOperands, std::vector<int> inCoeffs);
    virtual Symbol* evaluate(Symbol* vars[]);
    virtual std::stringstream print();
    virtual AdditionOp* diff(int varIndex);
    
    virtual void combineAdds();
    virtual void combineMults();
    virtual void cleanMults();
    virtual void cleanAdds();
    
    virtual void addOperand(Operand* op) {if(op != 0) {operands.push_back(op); numOperands++; coeffs.push_back(1);}};
    virtual void addOperand(Operand* op, int coeff) {if(op != 0) {operands.push_back(op);
                                                numOperands++; coeffs.push_back(coeff);}};
    
    std::vector<int> getCoeffs() {return coeffs;};
    
    
};



#endif /* defined(__Bertini2__AdditionOp__) */
