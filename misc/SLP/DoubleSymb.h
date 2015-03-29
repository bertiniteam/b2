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
#include "Symbol.h"


class DoubleSymb : Symbol
{
public:
    virtual Symbol* addSymb(Symbol* operand);
    virtual Symbol* multSymb(Symbol* operand);
    virtual Symbol* expSym(int exp);
   
};

#endif /* defined(__Bertini2__DoubleSymb__) */
