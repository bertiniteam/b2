//
//  Test_Diff.cpp
//  Bertini2
//
//  Created by James Collins on 5/1/14.
//  Copyright (c) 2014 James Collins. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <typeinfo>
#include "AdditionOp.h"
#include "MultOp.h"
#include "ExpOp.h"
#include "LeafOperand.h"
#include "Factor.h"
#include "DoubleSymb.h"
#include "IntSymb.h"

using namespace boost::numeric::ublas;

int main()
{
    Symbol *varVal[3];
    
    varVal[0] = new DoubleSymb(2.1);
    varVal[1] = new DoubleSymb(-1);
    varVal[2] = new DoubleSymb(3.7);

    LeafOperand* c2 = new LeafOperand(new IntSymb(2));
    
    LeafOperand* x = new LeafOperand(0);
    LeafOperand* y = new LeafOperand(1);
    LeafOperand* z = new LeafOperand(2);
    
    ExpOp* xe2 = new ExpOp(x,2);
    ExpOp* xe3 = new ExpOp(x,3);
    ExpOp* ye2 = new ExpOp(y,2);
    ExpOp* ye3 = new ExpOp(y,3);
    ExpOp* ze2 = new ExpOp(z,2);



    // p(x,y)
    AdditionOp* p = new AdditionOp();
    AdditionOp* inner = new AdditionOp();
    MultOp* temp = new MultOp();
    
    inner->addOperand(x);
    inner->addOperand(z);
    
    temp->addOperand(inner);
    temp->addOperand(xe3);
    p->addOperand(temp);
    
//    p->addOperand(ye2);

//    p->addOperand(ze2);






    std::cout << "**************p(x,y)*****************" << std::endl;
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
    std::cout << p->print().str() << std::endl;
    std::cout << "!!!!!!!Differentiate!!!!!!!!!!!" << std::endl;
    p = p->diff(2);
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
    std::cout << p->print().str() << std::endl;
    Operand::clean(p);
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
    std::cout << p->print().str() << std::endl;
    
    
    
    
    
    
    
    
    
    return 1;
    
    
}
