//
//  Operand.h
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__Operand__
#define __Bertini2__Operand__

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "Symbol.h"


/*******************Operand***************************
Desc: The super class for all nodes in the polynomial tree.  Contains 
 every method that is used recursively in the tree, such as evaluate,
 print and diff.
*****************************************************/






class Operand
{

protected:
    //List of operands.  Used only for Addition and Mult Ops
    std::vector<Operand*> operands;
    //Number of operands, i.e. size of operands vector
    int numOperands;
    
    //Type of operand, add, mult, exp, leaf
    int type;
    
public:
    //Constants representing various Types
    static const int TYPE_ADD = 1;
    static const int TYPE_MULT = 2;
    static const int TYPE_EXP = 3;
    static const int TYPE_LEAF = 4;
    
    
    
    /////////////////Abstract Methods///////////////////////////////
    
    //How we evaluate the polynomial
    virtual Symbol* evaluate(Symbol* vars[]) = 0;
    //output the polynomial for debugging
    virtual std::stringstream print() = 0;
    //Differentiate the polynomial.  Output is dirty.
    virtual Operand* diff(int varIndex) = 0;
    
    //The cleaning methods
    virtual void combineAdds() = 0;
    virtual void combineMults() = 0;
    virtual void cleanMults() = 0;
    virtual void cleanAdds() = 0;
    static void clean(Operand* node) {node->combineAdds(); node->combineMults(); node->cleanMults(); node->cleanAdds();};
    ////////////////////////////////////////////////////////////////
    
    
    /////////Get Methods/////////////////
    int getType() {return type;};
    virtual std::vector<Operand*> getOperands() {return operands;};
    virtual int getNumOperands() {return numOperands;};
    
    //Manipulate operands vector
    virtual void addOperand(Operand* op) {if(op != 0) {operands.push_back(op); numOperands++;}};
    void clear() {operands.clear(); numOperands = 0;};
    
    
};


#endif /* defined(__Bertini2__Operand__) */
