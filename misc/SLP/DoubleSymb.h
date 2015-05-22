//
//  DoubleSymb.h
//  Bertini2
//
//  Created by James Collins on 12/9/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__DoubleSymb__
#define __Bertini2__DoubleSymb__

#include <iostream>
#include <math.h>
#include "Symbol.h"
#include "IntSymb.h"


class DoubleSymb : public Symbol
{
private:
    double value;
    
public:
    
    
    DoubleSymb() {};
    DoubleSymb(DoubleSymb const & copy) : value(copy.value) {isVar = copy.isVar;};
    DoubleSymb(double inValue) {value = inValue; isVar = 0;};
    
    virtual std::stringstream print();
    virtual Symbol* add(Symbol* operand);
    virtual Symbol* sub(Symbol* operand);
    virtual Symbol* mult(Symbol* operand);
    virtual Symbol* exp(int exp);
    virtual DoubleSymb* neg();
    
    virtual DoubleSymb* clone ()
    {
        return new DoubleSymb(*this);
    }
   
    double getValue() {return value;};
};

#endif /* defined(__Bertini2__DoubleSymb__) */
