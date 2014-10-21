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
#include <sstream>
#include <typeinfo>


class Symbol
{
protected:
    bool isVar;
    
public:
    
    virtual std::stringstream print() = 0;
    
    virtual Symbol* add(Symbol* operand) = 0;
    virtual Symbol* sub(Symbol* operand) = 0;
    virtual Symbol* mult(Symbol* operand) = 0;
    virtual Symbol* exp(int exp) = 0;
    virtual Symbol* neg() = 0;
    
    virtual Symbol* clone() = 0;
    
    
    
    
    virtual ~Symbol() {};
};


#include "DoubleSymb.h"
#include "IntSymb.h"
//#include "amp.hpp"

#endif /* defined(__Bertini2__Symbol__) */
