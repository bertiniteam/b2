//
//  Symbol.h
//  Bertini2
//
//  Created by James Collins on 12/9/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__Symbol__
#define __Bertini2__Symbol__

#include <iostream>

class Symbol
{
private:
    bool isVar;
    
public:
    Symbol() {};
    virtual Symbol* addSymb(Symbol* operand) = 0;
    virtual Symbol* multSymb(Symbol* operand) = 0;
    virtual Symbol* expSym(int exp) = 0;
};

#endif /* defined(__Bertini2__Symbol__) */
