//
//  Test_TemplateLeaf.cpp
//  Bertini2
//
//  Created by James Collins on 4/10/14.
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
    double xd = -12.1;
    double yd = 1;
    double zd = 1;
    
    double pexact = (zd + xd * yd + 10) * (xd * zd + yd + 30) * (yd * zd + xd + 20);
    
    Symbol *varVal[3];
    
    varVal[0] = new DoubleSymb(xd);
    varVal[1] = new DoubleSymb(yd);
    varVal[2] = new DoubleSymb(zd);
    
    
    LeafOperand* c21 = new LeafOperand(new IntSymb(21));
    LeafOperand* c320 = new LeafOperand(new IntSymb(320));
    LeafOperand* c30 = new LeafOperand(new IntSymb(30));
    LeafOperand* c600 = new LeafOperand(new IntSymb(600));
    LeafOperand* c20 = new LeafOperand(new IntSymb(20));
    LeafOperand* c40 = new LeafOperand(new IntSymb(40));
    LeafOperand* c810 = new LeafOperand(new IntSymb(810));
    LeafOperand* c10 = new LeafOperand(new IntSymb(10));
    LeafOperand* c200 = new LeafOperand(new IntSymb(200));
    LeafOperand* c300 = new LeafOperand(new IntSymb(300));
    LeafOperand* c6000 = new LeafOperand(new IntSymb(6000));

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
    MultOp* temp = new MultOp();
    
    temp->addOperand(x);
    temp->addOperand(ye2);
    temp->addOperand(ze2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
//    temp->addOperand(c21);
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp,21);
   
    temp = new MultOp();
    temp->addOperand(ye2);
    temp->addOperand(ze2);
    p->addOperand(temp);

    temp = new MultOp();
//    temp->addOperand(c320);
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp,320);
    
    temp = new MultOp();
//    temp->addOperand(c30);
    temp->addOperand(y);
    temp->addOperand(ze2);
    p->addOperand(temp,30);
    
    temp = new MultOp();
//    temp->addOperand(c30);
    temp->addOperand(x);
    temp->addOperand(z);
    p->addOperand(temp,30);
    
    temp = new MultOp();
//    temp->addOperand(c600);
    temp->addOperand(z);
    p->addOperand(temp,600);
    
    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(ye3);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(ye2);
    p->addOperand(temp);
    
    temp = new MultOp();
//    temp->addOperand(c21);
    temp->addOperand(xe2);
    temp->addOperand(ye2);
    p->addOperand(temp,21);
    
    temp = new MultOp();
    temp->addOperand(x);
    temp->addOperand(ye3);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
//    temp->addOperand(c20);
    temp->addOperand(x);
    temp->addOperand(ye2);
    p->addOperand(temp,20);
    
    temp = new MultOp();
//    temp->addOperand(c40);
    temp->addOperand(x);
    temp->addOperand(ye2);
    temp->addOperand(z);
    p->addOperand(temp,40);
    
    temp = new MultOp();
//    temp->addOperand(c40);
    temp->addOperand(xe2);
    temp->addOperand(y);
    p->addOperand(temp,40);
    
    temp = new MultOp();
//    temp->addOperand(c810);
    temp->addOperand(x);
    temp->addOperand(y);
    p->addOperand(temp,810);
    
    temp = new MultOp();
//    temp->addOperand(c10);
    temp->addOperand(ye2);
    temp->addOperand(z);
    p->addOperand(temp,10);
    
    temp = new MultOp();
//    temp->addOperand(c200);
    temp->addOperand(y);
    p->addOperand(temp,200);
    
    temp = new MultOp();
//    temp->addOperand(c300);
    temp->addOperand(x);
    p->addOperand(temp,300);
    
    temp = new MultOp();
//    temp->addOperand(c6000);
    p->addOperand(temp,6000);















    
    




    
    int pre_mults = 0;
    int factor_mults = 0;
    double peval;
    
    
    std::cout << "************** Appendix 1 Test Poly *****************" << std::endl;
    std::cout << "************** Evaluate Poly *****************" << std::endl;
    peval = ((DoubleSymb*)p->evaluate(varVal))->getValue();
    std::cout << "Unfactored Error = " << abs(peval-pexact) << std::endl;
    std::cout << p->print().str() << std::endl;
    pre_mults += Factor::countMults(p);
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;
    
    std::cout << "!!!!Factoring!!!!" << std::endl;
    p = Factor::factorAdd(p);
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
    std::cout << p->print().str() << std::endl;
    factor_mults += Factor::countMults(p);
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;

    
    std::cout << "************** Test Differentiate *****************" << std::endl;
    Operand* p_x = p->diff(0);
    Symbol* e_x = p->evaluate(varVal)->sub(p_x->evaluate(varVal));
    std::cout << "|p-p_x| = " << e_x->print().str() << std::endl;
  
    
    return 1;


}