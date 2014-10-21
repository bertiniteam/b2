//
//  IntSymbol.h
//  Bertini2
//
//  Created by James Collins on 12/24/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__IntSymbol__
#define __Bertini2__IntSymbol__

#include <iostream>
#include <math.h>
#include "Symbol.h"


class IntSymb : public Symbol
{
private:
    int value;
    
public:
    
    
    IntSymb() {};
    IntSymb(IntSymb const & copy) : value(copy.value) {isVar = copy.isVar;};
    IntSymb(int inValue) {value = inValue; isVar = 0;};
    
    virtual std::stringstream print();
    virtual Symbol* add(Symbol* operand);
    virtual Symbol* sub(Symbol* operand);
    virtual Symbol* mult(Symbol* operand);
    virtual IntSymb* exp(int exp);
    virtual IntSymb* neg();
    
    IntSymb* clone ()
    {
        return new IntSymb(*this);
    }
    
    int getValue() {return value;};
};

#endif /* defined(__Bertini2__IntSymbol__) */
