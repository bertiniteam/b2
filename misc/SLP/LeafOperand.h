//
//  LeafOperand.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__LeafOperand__
#define __Bertini2__LeafOperand__

#include <iostream>
#include "Operand.h"
#include "Symbol.h"


class LeafOperand : public Operand
{
private:
    //Leaf specific variables
    bool isVar; //is it a variable
    Symbol* symb; //If it's a coefficient, this stores the value
    int varIndex; //if a variable, this is the index that specifies which variable
    char* varName;
    bool isZero, isOne;  //Allows us to specify Leaf as specific value of 1 or 0.
    
public:
    LeafOperand() {type = TYPE_LEAF;};
    LeafOperand(Symbol* inSymb) {symb = inSymb; type = TYPE_LEAF; isVar = 0;};
    LeafOperand(int index) {varIndex = index; type = TYPE_LEAF; isVar = 1;};
    LeafOperand(bool one, bool zero) {type = TYPE_LEAF; isVar = 0; isZero = zero; isOne = one;
        if(zero)
        {
            symb = new IntSymb(0);
        }
        else if(one)
        {
            symb = new IntSymb(1);
        }};
    
    
    virtual Symbol* evaluate(Symbol* vars[]);
    virtual std::stringstream print();
    virtual LeafOperand* diff(int invarIndex);
    
    virtual void combineAdds() {};
    virtual void combineMults() {};
    virtual void cleanMults() {};
    virtual void cleanAdds() {};
        
    bool getIsVar() {return isVar;};
    Symbol* getSymb() {return symb;};
    bool getIsZero() {return isZero;};
    bool getIsOne() {return isOne;};
    
    virtual void addOperand(Operand* op) {};
};

#endif /* defined(__Bertini2__LeafOperand__) */
