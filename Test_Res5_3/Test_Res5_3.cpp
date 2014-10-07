//
//  Test_Res5_3.cpp
//  Bertini2
//
//  Created by James Collins on 4/29/14.
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
    Symbol *varVal[10];
    
    varVal[0] = new DoubleSymb(1.2);
    varVal[1] = new DoubleSymb(-2.1);
    varVal[2] = new DoubleSymb(1.7);
    varVal[3] = new DoubleSymb(1.7);
    varVal[4] = new DoubleSymb(0.2);
    varVal[5] = new DoubleSymb(-1.1);
    varVal[6] = new DoubleSymb(1.7);
    varVal[7] = new DoubleSymb(1.7);
    varVal[8] = new DoubleSymb(-1.1);
    varVal[9] = new DoubleSymb(1.7);


    
    
    
    LeafOperand* c2 = new LeafOperand(new IntSymb(2));
    LeafOperand* c3 = new LeafOperand(new IntSymb(3));
    LeafOperand* c4 = new LeafOperand(new IntSymb(4));
    LeafOperand* c5 = new LeafOperand(new IntSymb(5));
    LeafOperand* c6 = new LeafOperand(new IntSymb(6));
    LeafOperand* c7 = new LeafOperand(new IntSymb(7));
    
    
    LeafOperand* a0 = new LeafOperand(0);
    LeafOperand* a1 = new LeafOperand(1);
    LeafOperand* a2 = new LeafOperand(2);
    LeafOperand* a3 = new LeafOperand(3);
    LeafOperand* a4 = new LeafOperand(4);
    LeafOperand* a5 = new LeafOperand(5);
    
    LeafOperand* b0 = new LeafOperand(6);
    LeafOperand* b1 = new LeafOperand(7);
    LeafOperand* b2 = new LeafOperand(8);
    LeafOperand* b3 = new LeafOperand(9);

    
    ExpOp *a0e[6], *a1e[6], *a2e[6], *a3e[6], *a4e[6], *a5e[6];
    ExpOp *b0e[6], *b1e[6], *b2e[6], *b3e[6];
    
    for(int ii = 2; ii <= 5; ii++)
    {
        a0e[ii] = new ExpOp(a0,ii);
        a1e[ii] = new ExpOp(a1,ii);
        a2e[ii] = new ExpOp(a2,ii);
        a3e[ii] = new ExpOp(a3,ii);
        a4e[ii] = new ExpOp(a4,ii);
        a5e[ii] = new ExpOp(a5,ii);
        
        b0e[ii] = new ExpOp(b0,ii);
        b1e[ii] = new ExpOp(b1,ii);
        b2e[ii] = new ExpOp(b2,ii);
        b3e[ii] = new ExpOp(b3,ii);
    }
    
    
    
    
    
    
    // p(x,y)
    AdditionOp* p = new AdditionOp();
    MultOp* temp = new MultOp();
    
    temp->addOperand(c5);
    temp->addOperand(b3);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0e[3]);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0e[3]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0e[4]);
    temp->addOperand(a4);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[3]);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0);
    temp->addOperand(a0e[2]);
    p->addOperand(temp,-1);

    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[3]);
    temp->addOperand(a3e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(a5e[2]);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[3]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[2]);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0);
    temp->addOperand(a0e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[3]);
    temp->addOperand(a3e[2]);
    p->addOperand(temp,1);
    
    //Term 10
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a4);
    temp->addOperand(b0);
    temp->addOperand(b2);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(b1e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0e[2]);
    temp->addOperand(b1);
    temp->addOperand(a3);
    temp->addOperand(b2);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3);
    temp->addOperand(a0e[2]);
    temp->addOperand(b1);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0);
    temp->addOperand(a2);
    temp->addOperand(a4);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0);
    temp->addOperand(a5);
    temp->addOperand(b0e[2]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3);
    temp->addOperand(a0);
    temp->addOperand(a5e[2]);
    temp->addOperand(b2);
    temp->addOperand(b0e[3]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(a5);
    temp->addOperand(b0e[3]);
    temp->addOperand(a4);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(a4e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    //Term 20
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[4]);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(a2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a4e[2]);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a0 );
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a4e[2]);
    temp->addOperand(b2);
    temp->addOperand(b0e[3]);
    temp->addOperand(a2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a4e[2]);
    temp->addOperand(b0e[4]);
    temp->addOperand(a5);
    temp->addOperand(b2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b0e[3]);
    temp->addOperand(b2);
    temp->addOperand(a3e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[2]);
    temp->addOperand(a3e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b0e[3]);
    temp->addOperand(a2);
    temp->addOperand(a4);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b0e[3]);
    temp->addOperand(a5);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[4]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a3);
    temp->addOperand(b0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[4]);
    temp->addOperand(a0e[2]);
    temp->addOperand(b1);
    temp->addOperand(a2);
    p->addOperand(temp,1);
    

    //Term 30
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a4);
    temp->addOperand(b1e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0);
    temp->addOperand(a3e[2]);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b1e[2]);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4e[2]);
    temp->addOperand(b1e[4]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a1e[2]);
    temp->addOperand(a4);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[4]);
    temp->addOperand(a1e[2]);
    temp->addOperand(a0);
    temp->addOperand(b1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[4]);
    temp->addOperand(a0e[2]);
    temp->addOperand(a1);
    temp->addOperand(b2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a4e[2]);
    temp->addOperand(a1);
    temp->addOperand(b0e[3]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3e[2]);
    temp->addOperand(a0);
    temp->addOperand(b1e[3]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[5]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    //Term 40
    temp = new MultOp();
    temp->addOperand(b2e[2]);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2);
    temp->addOperand(a0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b2e[3]);
    temp->addOperand(a0e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(a0e[2]);
    temp->addOperand(b2e[4]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[5]);
    temp->addOperand(a0e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0e[3]);
    temp->addOperand(a5);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(b0e[4]);
    temp->addOperand(a5e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[5]);
    temp->addOperand(a0e[3]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a1);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b1e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a3);
    temp->addOperand(b1e[3]);
    temp->addOperand(b0);
    temp->addOperand(a1);
    temp->addOperand(a5);
    p->addOperand(temp,-1);
    
    //Term 50
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b2e[4]);
    temp->addOperand(a2);
    temp->addOperand(a0);
    temp->addOperand(a5);
    temp->addOperand(b0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b2e[3]);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(a0);
    temp->addOperand(a4);
    temp->addOperand(b0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b2e[2]);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2);
    temp->addOperand(a0);
    temp->addOperand(a3);
    temp->addOperand(b0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2);
    temp->addOperand(b1);
    temp->addOperand(a1);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1);
    temp->addOperand(a1);
    temp->addOperand(a0);
    p->addOperand(temp,-1);

    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b2e[3]);
    temp->addOperand(b1);
    temp->addOperand(a1);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[4]);
    temp->addOperand(b1);
    temp->addOperand(a1);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a3e[4]);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a3);
    temp->addOperand(b0e[3]);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b0e[2]);
    temp->addOperand(b1);
    temp->addOperand(a4);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    //Term 60
    temp = new MultOp();
    temp->addOperand(b3e[3]);
    temp->addOperand(a3);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(a2);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(b1e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b1e[2]);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b1e[3]);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(a2);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(b1e[2]);
    temp->addOperand(a2);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b1e[2]);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b1);
    temp->addOperand(b3);
    temp->addOperand(a2e[2]);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b1e[2]);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0);
    temp->addOperand(a2);
    temp->addOperand(a4);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a4);
    temp->addOperand(b1e[4]);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    //Term 70
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b1e[3]);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[3]);
    temp->addOperand(a2);
    temp->addOperand(a4);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a2);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c7);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[3]);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(b1);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[3]);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(b2e[2]);
    temp->addOperand(a4);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    //Term 80
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(b2);
    temp->addOperand(a4e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(b1);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(a3);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3);
    temp->addOperand(a4e[2]);
    temp->addOperand(b2);
    temp->addOperand(b0);
    temp->addOperand(b1e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b0e[3]);
    temp->addOperand(a5);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(a2);
    temp->addOperand(b0e[2]);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(b2e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[3]);
    temp->addOperand(b2);
    temp->addOperand(b0);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0e[3]);
    temp->addOperand(b2);
    temp->addOperand(a2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a4);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    //Term 90
    temp = new MultOp();
    temp->addOperand(c7);
    temp->addOperand(b3e[2]);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(a0);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[2]);
    temp->addOperand(a5);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1e[3]);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[3]);
    temp->addOperand(a2);
    temp->addOperand(a4);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c5);
    temp->addOperand(b3e[2]);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(a2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a0);
    temp->addOperand(a4);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1e[3]);
    temp->addOperand(a3);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(a3);
    temp->addOperand(b0e[4]);
    temp->addOperand(a5e[2]);
    temp->addOperand(b2);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(a3);
    temp->addOperand(b0e[4]);
    temp->addOperand(a5);
    temp->addOperand(a4);
    p->addOperand(temp,1);
    
    //Term 100
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a3);
    temp->addOperand(b0e[3]);
    temp->addOperand(a4e[2]);
    temp->addOperand(b1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3e[3]);
    temp->addOperand(a3);
    temp->addOperand(b0e[2]);
    temp->addOperand(a2);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1e[3]);
    temp->addOperand(a2);
    temp->addOperand(a5e[2]);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b1e[3]);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2);
    temp->addOperand(a4);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b1);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2e[2]);
    temp->addOperand(a4);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1e[2]);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(a4e[2]);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2);
    temp->addOperand(a3e[2]);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4e[2]);
    temp->addOperand(b1e[3]);
    temp->addOperand(b0);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3e[2]);
    temp->addOperand(a1e[2]);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,1);
    
    //Term 110
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(b1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[2]);
    temp->addOperand(a1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b1e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a1);
    temp->addOperand(a3);
    temp->addOperand(b1e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[2]);
    temp->addOperand(a1);
    temp->addOperand(a5);
    temp->addOperand(b1e[3]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3e[2]);
    temp->addOperand(a3e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b3);
    temp->addOperand(a3);
    temp->addOperand(b1e[4]);
    temp->addOperand(a5);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5e[2]);
    temp->addOperand(b1e[4]);
    temp->addOperand(b0);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b2e[3]);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a5);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b2e[2]);
    temp->addOperand(b3);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a4);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    //Term 120
    temp = new MultOp();
    temp->addOperand(c2);
    temp->addOperand(b2);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b2);
    temp->addOperand(b3e[3]);
    temp->addOperand(a2);
    temp->addOperand(b0);
    temp->addOperand(a1e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0);
    temp->addOperand(a1e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b2e[3]);
    temp->addOperand(b0);
    temp->addOperand(a1e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[4]);
    temp->addOperand(b0);
    temp->addOperand(a1e[2]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[2]);
    temp->addOperand(a3e[3]);
    temp->addOperand(b0e[3]);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3e[4]);
    temp->addOperand(a1e[3]);
    temp->addOperand(b0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4e[3]);
    temp->addOperand(b0e[4]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b3e[3]);
    temp->addOperand(a2e[3]);
    temp->addOperand(b0e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a5e[3]);
    temp->addOperand(b0e[5]);
    p->addOperand(temp,1);
    
    //Term 130
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(a1);
    temp->addOperand(a0);
    p->addOperand(temp,1);

    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1);
    temp->addOperand(b0);
    temp->addOperand(a1e[2]);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c6);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(a2);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c7);
    temp->addOperand(b3e[2]);
    temp->addOperand(a0);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(b1);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c7);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(a0);
    temp->addOperand(a4);
    temp->addOperand(b2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[3]);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(a3);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c6);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(a5);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(a5);
    temp->addOperand(b1);
    temp->addOperand(b0e[2]);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    //Term 140
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(a3);
    temp->addOperand(a0);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(b3e[2]);
    temp->addOperand(a2);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(a3);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(c4);
    temp->addOperand(b3e[2]);
    temp->addOperand(a4);
    temp->addOperand(b2);
    temp->addOperand(b0);
    temp->addOperand(b1);
    temp->addOperand(a2);
    temp->addOperand(a0);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b0e[2]);
    temp->addOperand(a3);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(b1e[2]);
    temp->addOperand(a2);
    temp->addOperand(a5);
    temp->addOperand(b0e[2]);
    temp->addOperand(a4);
    temp->addOperand(b2);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(c3);
    temp->addOperand(b3);
    temp->addOperand(b1e[2]);
    temp->addOperand(a2);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b3);
    temp->addOperand(a4);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    temp = new MultOp();
    temp->addOperand(a4);
    temp->addOperand(b1e[3]);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b2);
    temp->addOperand(a1);
    p->addOperand(temp,-1);
    
    //Term 150
    temp = new MultOp();
    temp->addOperand(a5);
    temp->addOperand(b2e[2]);
    temp->addOperand(b1e[2]);
    temp->addOperand(b0);
    temp->addOperand(a3);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(b3);
    temp->addOperand(a2);
    temp->addOperand(b0);
    temp->addOperand(a4);
    temp->addOperand(b2e[2]);
    temp->addOperand(a1);
    p->addOperand(temp,1);
    
    temp = new MultOp();
    temp->addOperand(b1);
    temp->addOperand(a2);
    temp->addOperand(b0);
    temp->addOperand(a5);
    temp->addOperand(b2e[3]);
    temp->addOperand(a1);
    p->addOperand(temp,-1);


    
    
    
    
    
    
    int pre_mults = 0;
    int factor_mults = 0;
    
    
    std::cout << "**************p(x,y)*****************" << std::endl;
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
//    std::cout << p->print().str() << std::endl;
    pre_mults += Factor::countMults(p);
    std::cout << "+'s = " << Factor::countAdds(p) << std::endl;
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;
    
    std::cout << "!!!!Factoring!!!!" << std::endl;
    p = Factor::factorAdd(p);
    std::cout << "Evaluate = " << p->evaluate(varVal)->print().str() << std::endl;
//    std::cout << p->print().str() << std::endl;
    factor_mults += Factor::countMults(p);
    std::cout << "*'s = " << Factor::countMults(p) << std::endl;
    
    
    
    return 1;
    
}