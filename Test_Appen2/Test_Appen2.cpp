//
//  Test_Appen2.cpp
//  Bertini2
//
//  Created by James Collins on 4/25/14.
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
    double xd = 4.34;
    double yd = 4.1;
    double zd = -2.333;
    
    double pexact = (pow(xd, 0.3e1) * (zd + yd) + zd - 11) * (xd * xd * (zd * zd + yd * yd) + yd + 90);
    double pxexact = 3 * xd * xd * (zd + yd) * (xd * xd * (zd * zd + yd * yd) + yd + 90) + 2 * ( pow((double) xd, (double) 3) * (zd + yd) + zd - 11) * xd * (zd * zd + yd * yd);
    double pyexact = pow(xd, 0.3e1) * (xd * xd * (zd * zd + yd * yd) + yd + 0.90e2) + (pow(xd, 0.3e1) * (zd + yd) + zd - 0.11e2) * (0.2e1 * xd * xd * yd + 0.1e1);
    double pzexact = (pow(xd, 0.3e1) + 0.1e1) * (xd * xd * (zd * zd + yd * yd) + yd + 0.90e2) + 0.2e1 * (pow(xd, 0.3e1) * (zd + yd) + zd - 0.11e2) * xd * xd * zd;
    
    
    
    
    Symbol *varVal[3];
    
    varVal[0] = new DoubleSymb(xd);
    varVal[1] = new DoubleSymb(yd);
    varVal[2] = new DoubleSymb(zd);
    
    
    LeafOperand* c90 = new LeafOperand(new IntSymb(90));
    LeafOperand* c11 = new LeafOperand(new IntSymb(11));
    LeafOperand* c990 = new LeafOperand(new IntSymb(990));
    
    LeafOperand* x = new LeafOperand(0);
    LeafOperand* y = new LeafOperand(1);
    LeafOperand* z = new LeafOperand(2);
    
    ExpOp* xe2 = new ExpOp(x,2);
    ExpOp* xe3 = new ExpOp(x,3);
    ExpOp* xe4 = new ExpOp(x,4);
    ExpOp* xe5 = new ExpOp(x,5);
    ExpOp* ye2 = new ExpOp(y,2);
    ExpOp* ye3 = new ExpOp(y,3);
    ExpOp* ye4 = new ExpOp(y,4);
    ExpOp* ye5 = new ExpOp(y,5);
    ExpOp* ze2 = new ExpOp(z,2);
    ExpOp* ze3 = new ExpOp(z,3);
    ExpOp* ze4 = new ExpOp(z,4);
    ExpOp* ze5 = new ExpOp(z,5);
    
    
    // p(x,y)
    AdditionOp* p = new AdditionOp();
    MultOp* temp = new MultOp();
    
    temp->addOperand(xe5);
    temp->addOperand(ze3);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe5);
    temp->addOperand(ye2);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(c90);
    temp->addOperand(xe3);
    temp->addOperand(z);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(xe5);
    temp->addOperand(y);
    temp->addOperand(ze2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe5);
    temp->addOperand(ye3);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(ye2);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(c90);
    temp->addOperand(xe3);
    temp->addOperand(y);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(ze3);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(ye2);
    temp->addOperand(z);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(c90);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(c11);
    temp->addOperand(xe2);
    temp->addOperand(ze2);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c11);
    temp->addOperand(xe2);
    temp->addOperand(ye2);
    p->addOperand(temp,-1);

    temp = new MultOp();
    temp->addOperand(c11);
    temp->addOperand(y);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c990);
    p->addOperand(temp,-1);



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    int pre_mults = 0;
    int factor_mults = 0;
    double peval, pxeval, pyeval, pzeval;
    
    
    std::cout << "************** Appendix 2 Test Poly *****************" << std::endl;
    std::cout << "************** Unfactored Tests *****************" << std::endl;
    peval = ((DoubleSymb*)p->evaluate(varVal))->getValue();
    std::cout << "Unfactored Error = " << abs(pexact-peval) << std::endl;
    std::cout << p->print().str() << std::endl;
    pre_mults += Factor::countMults(p);
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;
    Operand* p_x = p->diff(0);
    pxeval = ((DoubleSymb*)p_x->evaluate(varVal))->getValue();
    std::cout << "Diff_x Error = " << abs(pxeval-pxexact) << std::endl;
    Operand* p_y = p->diff(1);
    pyeval = ((DoubleSymb*)p_y->evaluate(varVal))->getValue();
    std::cout << "Diff_y Error = " << abs(pyeval-pyexact) << std::endl;
    Operand* p_z = p->diff(2);
    pzeval = ((DoubleSymb*)p_z->evaluate(varVal))->getValue();
    std::cout << "Diff_z Error = " << abs(pzeval-pzexact) << std::endl;
    
    std::cout << "************** Factored Tests *****************" << std::endl;
    std::cout << "!!!!!!!!!!!!Factoring!!!!!!!!!!!!!!" << std::endl;
    p = Factor::factorAdd(p);
    peval = ((DoubleSymb*)p->evaluate(varVal))->getValue();
    std::cout << "Factored Error = " << abs(peval-pexact) << std::endl;
    std::cout << p->print().str() << std::endl;
    factor_mults += Factor::countMults(p);
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;
    
    p_x = p->diff(0);
    pxeval = ((DoubleSymb*)p_x->evaluate(varVal))->getValue();
    std::cout << "Diff_x Error = " << abs(pxeval-pxexact) << std::endl;
    p_y = p->diff(1);
    pyeval = ((DoubleSymb*)p_y->evaluate(varVal))->getValue();
    std::cout << "Diff_y Error = " << abs(pyeval-pyexact) << std::endl;
    p_z = p->diff(2);
    pzeval = ((DoubleSymb*)p_z->evaluate(varVal))->getValue();
    std::cout << "Diff_z Error = " << abs(pzeval-pzexact) << std::endl;

    
    
    
    return 1;
    
    
}