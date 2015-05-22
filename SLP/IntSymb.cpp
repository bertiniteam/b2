//
//  IntSymb.cpp
//  Bertini2
//
//  Created by James Collins on 12/24/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "IntSymb.h"


Symbol* IntSymb::add(Symbol* operand)
{
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        DoubleSymb* retop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb(tempop->getValue() + (double)this->value);
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        IntSymb* retop = new IntSymb(tempop->getValue() + this->value);
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





Symbol* IntSymb::sub(Symbol* operand)
{
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        DoubleSymb* retop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb((double)this->value - tempop->getValue());
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        IntSymb* retop = new IntSymb(this->value - tempop->getValue());
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









Symbol* IntSymb::mult(Symbol* operand)
{
    if(typeid(*operand) == typeid(DoubleSymb))
    {
        DoubleSymb* tempop;
        DoubleSymb* retop;
        tempop = (DoubleSymb*)operand;
        retop = new DoubleSymb(tempop->getValue() * (double)this->value);
        return retop;
    }
    else if(typeid(*operand) == typeid(IntSymb))
    {
        IntSymb* tempop = (IntSymb*)operand;
        IntSymb* retop = new IntSymb(tempop->getValue() * this->value);
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



IntSymb* IntSymb::exp(int exp)
{
    IntSymb* retop = new IntSymb((int)pow(this->value,exp));
    
    return retop;
    
}



IntSymb* IntSymb::neg()
{
    IntSymb* retop = new IntSymb(-(this->value));
    
    return retop;
    
}




std::stringstream IntSymb::print()
{
    std::stringstream ret;
    
    
    ret << value;
    
    return ret;
}
