//
//  DoubleSymb.cpp
//  Bertini2
//
//  Created by James Collins on 12/9/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "DoubleSymb.h"







Symbol* DoubleSymb::add(Symbol* operand)
{
    
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* retop, *tempop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb(tempop->getValue() + this->value);
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        DoubleSymb* retop = new DoubleSymb((double)tempop->getValue() + this->value);
        return retop;
    }
//    else if(typeid(*operand) == typeid(AMP))
//    {
//        AMP tempop = *((AMP*)operand);
//        AMP* retop = new AMP(this->value + tempop);
//        return retop;
//    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
}




Symbol* DoubleSymb::sub(Symbol* operand)
{
    
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* retop, *tempop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb(this->value - tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        DoubleSymb* retop = new DoubleSymb(this->value - (double)tempop->getValue());
        return retop;
    }
//    else if(typeid(*operand) == typeid(AMP))
//    {
//        AMP tempop = *((AMP*)operand);
//        AMP* retop = new AMP(this->value - tempop);
//        return retop;
//    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
}






Symbol* DoubleSymb::mult(Symbol* operand)
{
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* retop, *tempop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb(tempop->getValue() * this->value);
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        DoubleSymb* retop = new DoubleSymb((double)tempop->getValue() * this->value);
        return retop;
    }
//    else if(typeid(*operand) == typeid(AMP))
//    {
//        AMP tempop = *((AMP*)operand);
//        AMP* retop = new AMP(this->value * tempop);
//        return retop;
//    }
    else
    {
        Symbol* retop = 0;
        return retop;
    }
   
}





Symbol* DoubleSymb::exp(int exp)
{
    this->value = (double)pow(this->value,exp);
    
    return this;
    
}


DoubleSymb* DoubleSymb::neg()
{
    DoubleSymb* retop = new DoubleSymb(-(this->value));
    
    return retop;
    
}





std::stringstream DoubleSymb::print()
{
    std::stringstream ret;
    
    
    ret << value;
    
    return ret;
}












