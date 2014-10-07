//
//  Test_Appen13.cpp
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
    double xd = 1;
    double yd = 1;
    double zd = 1;
    double wd = 0;
    
    double pexact = (0.3e1 * pow(zd, 0.3e1) + 0.2e1 * wd * zd - 0.9e1 * pow(yd, 0.3e1) - yd * yd + 0.45e2 * pow(xd, 0.3e1)) * (wd * wd * pow(zd, 0.3e1) + 0.47e2 * xd * yd - wd * wd);
    double pxexact = 0.135e3 * (double) (xd * xd) * (wd * wd * pow(zd, 0.3e1) + (double) (47 * xd * yd) - wd * wd) + 0.47e2 * (0.3e1 * pow(zd, 0.3e1) + 0.2e1 * wd * zd - (double) (9 *  pow((double) yd, (double) 3)) - (double) (yd * yd) + (double) (45 *  pow((double) xd, (double) 3))) * (double) yd;

    double pyexact = (double) (-27 * yd * yd - 2 * yd) * (wd * wd * pow(zd, 0.3e1) + (double) (47 * xd * yd) - wd * wd) + 0.47e2 * (0.3e1 * pow(zd, 0.3e1) + 0.2e1 * wd * zd - (double) (9 *  pow((double) yd, (double) 3)) - (double) (yd * yd) + (double) (45 *  pow((double) xd, (double) 3))) * (double) xd;

    
    double pzexact = (9 * zd * zd + 2 * wd) * (wd * wd *  pow((double) zd, (double) 3) + 47 * xd * yd - wd * wd) + 3 * (3 *  pow((double) zd, (double) 3) + 2 * wd * zd - 9 *  pow((double) yd, (double) 3) - yd * yd + 45 *  pow((double) xd, (double) 3)) * wd * wd * zd * zd;

    double pwexact = 0.2e1 * zd * (wd * wd * pow(zd, 0.3e1) + (double) (47 * xd * yd) - wd * wd) + (0.3e1 * pow(zd, 0.3e1) + 0.2e1 * wd * zd - (double) (9.0 *  pow((double) yd, (double) 3)) - (double) (yd * yd) + (double) (45 *  pow((double) xd, (double) 3))) * (0.2e1 * wd * pow(zd, 0.3e1) - 0.2e1 * wd);


    
    Symbol *varVal[4];
    
    varVal[0] = new DoubleSymb(xd);
    varVal[1] = new DoubleSymb(yd);
    varVal[2] = new DoubleSymb(zd);
    varVal[3] = new DoubleSymb(wd);
    
    
    
    LeafOperand* c3 = new LeafOperand(new IntSymb(3));
    LeafOperand* c141 = new LeafOperand(new IntSymb(141));
    LeafOperand* cm3 = new LeafOperand(new IntSymb(-3));
    LeafOperand* c2 = new LeafOperand(new IntSymb(2));
    LeafOperand* c94 = new LeafOperand(new IntSymb(94));
    LeafOperand* cm2 = new LeafOperand(new IntSymb(-2));
    LeafOperand* cm9 = new LeafOperand(new IntSymb(-9));
    LeafOperand* cm423 = new LeafOperand(new IntSymb(-423));
    LeafOperand* c9 = new LeafOperand(new IntSymb(9));
    LeafOperand* cm1 = new LeafOperand(new IntSymb(-1));
    LeafOperand* cm47 = new LeafOperand(new IntSymb(-47));
    LeafOperand* c45 = new LeafOperand(new IntSymb(45));
    LeafOperand* c2115 = new LeafOperand(new IntSymb(2115));
    LeafOperand* cm45 = new LeafOperand(new IntSymb(-45));
    
    
    LeafOperand* x = new LeafOperand(0);
    LeafOperand* y = new LeafOperand(1);
    LeafOperand* z = new LeafOperand(2);
    LeafOperand* w = new LeafOperand(3);
    
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
    
    ExpOp* we2 = new ExpOp(w,2);
    ExpOp* we3 = new ExpOp(w,3);
    
    
    
    
    
    // p(x,y)
    AdditionOp* p = new AdditionOp();
    MultOp* temp = new MultOp();
    
    temp->addOperand(c3);
    temp->addOperand(we2);
    temp->addOperand(ze6);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(c141);
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(ze3);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm3);
    temp->addOperand(ze3);
    temp->addOperand(we2);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(ze4);
    temp->addOperand(we3);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(c94);
    temp->addOperand(x);
    temp->addOperand(y);
    temp->addOperand(z);
    temp->addOperand(w);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm2);
    temp->addOperand(z);
    temp->addOperand(we3);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm9);
    temp->addOperand(ye3);
    temp->addOperand(ze3);
    temp->addOperand(we2);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm423);
    temp->addOperand(x);
    temp->addOperand(ye4);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(c9);
    temp->addOperand(ye3);
    temp->addOperand(we2);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm1);
    temp->addOperand(ye2);
    temp->addOperand(ze3);
    temp->addOperand(we2);
    p->addOperand(temp);//
    
    temp = new MultOp();
    temp->addOperand(cm47);
    temp->addOperand(x);
    temp->addOperand(ye3);
    p->addOperand(temp);//
    
    temp = new MultOp();
//    temp->addOperand(c3);
    temp->addOperand(ye2);
    temp->addOperand(we2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(c45);
    temp->addOperand(xe3);
    temp->addOperand(ze3);
    temp->addOperand(we2);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(c2115);
    temp->addOperand(xe4);
    temp->addOperand(y);
    p->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(cm45);
    temp->addOperand(xe3);
    temp->addOperand(we2);
    p->addOperand(temp);













    
    
    
    
    int pre_mults = 0;
    int factor_mults = 0;
    double peval, pxeval, pyeval, pzeval, pweval;
    
    
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
    Operand* p_w = p->diff(3);
    pweval = ((DoubleSymb*)p_w->evaluate(varVal))->getValue();
    std::cout << "Diff_w Error = " << abs(pweval-pwexact) << std::endl;

    
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
    p_w = p->diff(3);
    pweval = ((DoubleSymb*)p_z->evaluate(varVal))->getValue();
    std::cout << "Diff_w Error = " << abs(pweval-pwexact) << std::endl;

//    std::cout << "************** Other Tests *****************" << std::endl;
//    std::cout << p_w->print().str() << std::endl;
   
    
    
    return 1;
    
}