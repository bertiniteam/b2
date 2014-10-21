#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
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
    std::vector<Operand*> store;
    std::vector<Operand*> terms;
    
    
    // Set values for x,y and z
    Symbol *x0,*x1,*x2, *x3, *x4, *x5, *x6, *x7;
    
    x0 = new DoubleSymb(1.3);
    x1 = new DoubleSymb(-5.2);
    x2 = new DoubleSymb(1);
    x3 = new DoubleSymb(3.1);
    x4 = new DoubleSymb(2.9);
    x5 = new DoubleSymb(1.3);
    x6 = new DoubleSymb(2);
    x7 = new DoubleSymb(4);
    
    Symbol* var[] = {x0,x1,x2,x3,x4,x5,x6};
    
    // Define coefficients
    Symbol *c1, *c2, *c3;
    
    c1 = new DoubleSymb(4.0);
    c2 = new DoubleSymb(6.0);
    c3 = new DoubleSymb(2.45);
    
    
    //Create leaves and suboperands for
    // f = c1((x_0)^2) + c2(x_1) + c3(x_0)(x_2)
//    LeafOperand* c1Coeff = new LeafOperand(c1);
//    LeafOperand* c2Coeff = new LeafOperand(c2);
//    LeafOperand* c3Coeff = new LeafOperand(c3);
//    
//    LeafOperand* x0Var = new LeafOperand(0);
//    LeafOperand* x1Var = new LeafOperand(1);
//    LeafOperand* x2Var = new LeafOperand(2);
//    ExpOp* x0e2 = new ExpOp(x0Var,2);
//    
//    store.push_back(c1Coeff);
//    store.push_back(x0e2);
//    terms.push_back(new MultOp(store));
//    store.erase(store.begin(),store.end());
//    
//    store.push_back(c2Coeff);
//    store.push_back(x1Var);
//    terms.push_back(new MultOp(store));
//    store.erase(store.begin(),store.end());
//
//    store.push_back(c3Coeff);
//    store.push_back(x0Var);
//    store.push_back(x2Var);
//    terms.push_back(new MultOp(store));
//    store.erase(store.begin(),store.end());
//
//    
//    AdditionOp* tree = new AdditionOp(terms);
 
    
    LeafOperand* c1Coeff = new LeafOperand(c1);
    LeafOperand* c2Coeff = new LeafOperand(c2);
    LeafOperand* c3Coeff = new LeafOperand(c3);

    LeafOperand* x0Var = new LeafOperand(0);
    LeafOperand* x1Var = new LeafOperand(1);
    LeafOperand* x2Var = new LeafOperand(2);
    
    LeafOperand* zVar = new LeafOperand(0);
    LeafOperand* uVar = new LeafOperand(1);
    LeafOperand* vVar = new LeafOperand(2);
    LeafOperand* aVar = new LeafOperand(3);
    LeafOperand* bVar = new LeafOperand(4);
    LeafOperand* qVar = new LeafOperand(5);
    LeafOperand* wVar = new LeafOperand(6);

    ExpOp* x0e2 = new ExpOp(x0Var,2);
    ExpOp* x0e3 = new ExpOp(x0Var,3);
    ExpOp* x0e4 = new ExpOp(x0Var,4);
    ExpOp* x0e5 = new ExpOp(x0Var,5);
    ExpOp* x0e6 = new ExpOp(x0Var,6);
    ExpOp* x0e7 = new ExpOp(x0Var,7);
    ExpOp* x1e2 = new ExpOp(x1Var,2);
    ExpOp* x1e3 = new ExpOp(x1Var,3);
    ExpOp* x1e4 = new ExpOp(x1Var,4);
    ExpOp* x1e5 = new ExpOp(x1Var,5);
    ExpOp* x1e6 = new ExpOp(x1Var,6);
    ExpOp* x1e7 = new ExpOp(x1Var,7);
    ExpOp* x2e2 = new ExpOp(x2Var,2);
    ExpOp* x2e3 = new ExpOp(x2Var,3);
    ExpOp* x2e4 = new ExpOp(x2Var,4);
    
    ExpOp* ue2 = new ExpOp(uVar,2);
    ExpOp* ue3 = new ExpOp(uVar,3);
    ExpOp* ue4 = new ExpOp(uVar,4);
    ExpOp* ve2 = new ExpOp(vVar,2);
    ExpOp* ve3 = new ExpOp(vVar,3);
    ExpOp* ve4 = new ExpOp(vVar,4);
    
    AdditionOp* x0px1 = new AdditionOp();
    x0px1->addOperand(x0Var);
    x0px1->addOperand(x1Var);
    
    
    
    
    
    AdditionOp* ex1 = new AdditionOp();
    
    MultOp* temp = new MultOp();
    temp->addOperand(zVar);
    ex1->addOperand(temp);
    

    
    
    
    AdditionOp* ex2 = new AdditionOp();
    temp = new MultOp();
    temp->addOperand(x0e2);
//    temp->addOperand(x1Var);
    ex2->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x0Var);
//    temp->addOperand(x1Var);
    ex2->addOperand(temp,2);
  
//    temp = new MultOp();
//    temp->addOperand(x0e2);
//    temp->addOperand(x1Var);
//    ex2->addOperand(temp);

    ex2->addOperand(x1Var,-1);

    
    
    ////////// Sum terms doesn't work for this example////////////
    AdditionOp* ex3 = new AdditionOp();
//    temp = new MultOp();
//    temp->addOperand(x0e7);
////    temp->addOperand(x1e7);
//    ex3->addOperand(temp);
//    
//    temp = new MultOp();
//    temp->addOperand(x0e6);
////    temp->addOperand(x1e7);
//    ex3->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x0e5);
//    temp->addOperand(x1e2);
    ex3->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x0e4);
//    temp->addOperand(x1e7);
    ex3->addOperand(temp);
    
    temp = new MultOp();
    temp->addOperand(x0e3);
    temp->addOperand(x1e5);
    ex3->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(x0e2);
    temp->addOperand(x1e6);
    ex3->addOperand(temp);

    temp = new MultOp();
    temp->addOperand(x0Var);
    temp->addOperand(x1e7);
    ex3->addOperand(temp);



    
    
    
    
    
    
    

    
    
    //Example 1
    std::cout << "Evaluate = " << ex1->evaluate(var)->print().str() << std::endl;
    std::cout << ex1->print().str() << std::endl;
    std::cout << "+'s = " << Factor::countAdds(ex1) << std::endl;
    std::cout << "*'s = " << Factor::countMults(ex1) << std::endl;

    ex1 = Factor::factorAdd(ex1);
    std::cout << "Evaluate = " << ex1->evaluate(var)->print().str() << std::endl;
    std::cout << ex1->print().str() << std::endl;
    std::cout << "+'s = " << Factor::countAdds(ex1) << std::endl;
    std::cout << "*'s = " << Factor::countMults(ex1) << std::endl;

    
    
    //Example 2
    std::cout << "Evaluate = " << ex2->evaluate(var)->print().str() << std::endl;
    std::cout << ex2->print().str() << std::endl;
    std::cout << "+'s = " << Factor::countAdds(ex2) << std::endl;
    std::cout << "*'s = " << Factor::countMults(ex2) << std::endl;
    
    ex2 = Factor::factorAdd(ex2);
    std::cout << "Evaluate = " << ex2->evaluate(var)->print().str() << std::endl;
    std::cout << ex2->print().str() << std::endl;
    std::cout << "+'s = " << Factor::countAdds(ex2) << std::endl;
    std::cout << "*'s = " << Factor::countMults(ex2) << std::endl;
    
    
    //Example 3
//    std::cout << "Evaluate = " << ex3->evaluate(var)->print().str() << std::endl;
//    std::cout << ex3->print().str() << std::endl;
//    std::cout << "+'s = " << Factor::countAdds(ex3) << std::endl;
//    std::cout << "*'s = " << Factor::countMults(ex3) << std::endl;
//    
//    ex3 = Factor::factorAdd(ex3);
//    std::cout << "Evaluate = " << ex3->evaluate(var)->print().str() << std::endl;
//    std::cout << ex3->print().str() << std::endl;
//    std::cout << "+'s = " << Factor::countAdds(ex3) << std::endl;
//    std::cout << "*'s = " << Factor::countMults(ex3) << std::endl;
    
    
    
    
    
    
    
    
  
    
    return 1;
}