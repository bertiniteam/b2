//
//  MultOp.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__MultOp__
#define __Bertini2__MultOp__

#include <iostream>
#include <vector>
#include "Operand.h"
#include "AdditionOp.h"


class MultOp : public Operand
{
public:
    MultOp() {operands = std::vector<Operand*>(); numOperands = 0; type = TYPE_MULT;};
    MultOp(std::vector<Operand*> inOperands);
    
    virtual Symbol* evaluate(Symbol* vars[]);
    virtual std::stringstream print();
    virtual AdditionOp* diff(int varIndex);
    
    virtual void combineAdds();
    virtual void combineMults();
    virtual void cleanMults();
    virtual void cleanAdds();
    
};

#endif /* defined(__Bertini2__MultOp__) */
