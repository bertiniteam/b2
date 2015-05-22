//
//  ExpOp.h
//  Bertini2
//
//  Created by James Collins on 11/20/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__ExpOp__
#define __Bertini2__ExpOp__

#include <iostream>
#include <vector>
#include "Operand.h"
#include "MultOp.h"


class ExpOp : public Operand
{
private:
    Operand* base;
    int exp;
public:
    ExpOp() {type = TYPE_EXP;};
    ExpOp(Operand* inBase, int inExp);
    
    virtual Symbol* evaluate(Symbol* vars[]);
    virtual std::stringstream print();
    virtual Operand* diff(int varIndex);
    
    virtual void combineAdds() {base->combineAdds();};
    virtual void combineMults() {base->combineMults();};
    virtual void cleanMults();
    virtual void cleanAdds() {base->cleanAdds();};
    
    Operand* getBase() {return base;};
    int getExp() {return exp;};
    void setExp(int inExp) {exp = inExp;};
    
    virtual void addOperand(Operand* op) {};
};

#endif /* defined(__Bertini2__ExpOp__) */
