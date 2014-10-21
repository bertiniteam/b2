//
//  Test_Appen3.cpp
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
    double xd = -3.8792878;
    double yd = 1.298877846;
    double zd = -12.899002792;
    
    double pexact = (yd * pow(zd, 0.3e1) + xd * yd * zd + yd * yd + pow(xd, 0.3e1)) * (xd * (pow(zd, 0.4e1) + 0.1e1) + zd + pow(xd, 0.3e1) * yd * yd);
    double pxexact = (yd * zd + 3 * xd * xd) * (xd * ( pow((double) zd, (double) 4) + 1) + zd +  pow((double) xd, (double) 3) * yd * yd) + (yd *  pow((double) zd, (double) 3) + xd * yd * zd + yd * yd +  pow((double) xd, (double) 3)) * ( pow((double) zd, (double) 4) + 1 + 3 * xd * xd * yd * yd);
    double pyexact = (pow(zd, 0.3e1) + xd * zd + (double) (2 * yd)) * (xd * (pow(zd, 0.4e1) + 0.1e1) + zd + pow(xd, 0.3e1) * (double) (yd * yd)) + 0.2e1 * ((double) yd * pow(zd, 0.3e1) + xd * (double) yd * zd + (double) (yd * yd) + pow(xd, 0.3e1)) * pow(xd, 0.3e1) * (double) yd;

    double pzexact = (3 * yd * zd * zd + xd * yd) * (xd * ( pow((double) zd, (double) 4) + 1) + zd +  pow((double) xd, (double) 3) * yd * yd) + (yd *  pow((double) zd, (double) 3) + xd * yd * zd + yd * yd +  pow((double) xd, (double) 3)) * (4 * xd *  pow((double) zd, (double) 3) + 1);

    
    
    
    
    Symbol *varVal[3];
    
    varVal[0] = new DoubleSymb(xd);
    varVal[1] = new DoubleSymb(yd);
    varVal[2] = new DoubleSymb(zd);
    
    
    
    LeafOperand* x = new LeafOperand(0);
    LeafOperand* y = new LeafOperand(1);
    LeafOperand* z = new LeafOperand(2);
    
    ExpOp* xe2 = new ExpOp(x,2);
    ExpOp* xe3 = new ExpOp(x,3);
    ExpOp* xe4 = new ExpOp(x,4);
    ExpOp* xe5 = new ExpOp(x,5);
    ExpOp* xe6 = new ExpOp(x,6);
    ExpOp* ye2 = new ExpOp(y,2);
    ExpOp* ye3 = new ExpOp(y,3);
    ExpOp* ye4 = new ExpOp(y,4);
    ExpOp* ye5 = new ExpOp(y,5);
    ExpOp* ze2 = new ExpOp(z,2);
    ExpOp* ze3 = new ExpOp(z,3);
    ExpOp* ze4 = new ExpOp(z,4);
    ExpOp* ze5 = new ExpOp(z,5);
    ExpOp* ze6 = new ExpOp(z,6);
    ExpOp* ze7 = new ExpOp(z,7);
    
    
    // p(x,y)
    AdditionOp* p = new AdditionOp();
    MultOp* temp = new MultOp();
    
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(ze7);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(ze3);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(y);
    temp->addOperand(ze4);
    p->addOperand(temp);
   
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(ye3);
    temp->addOperand(ze3);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(y);
    temp->addOperand(ze5);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe2);
    temp->addOperand(y);
    temp->addOperand(z);
    p->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(ze2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe4);
    temp->addOperand(ye3);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x);
    temp->addOperand(ye2);
    temp->addOperand(ze4);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x);
    temp->addOperand(ye2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(ye2);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(ye4);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe4);
    temp->addOperand(ze4);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe4);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe3);
    temp->addOperand(z);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(xe6);
    temp->addOperand(ye2);
    p->addOperand(temp);










    
    
    
    
    
    
    
    
    
    
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