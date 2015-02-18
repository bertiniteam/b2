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
#include "TreeIndex.h"



int main()
{

    double x=0,y=0,z=0, c1=1, c2=2, c3=3;
//    std::cout << "Enter x: ";
//    std::cin >> x;
//    std::cout << "Enter y: ";
//    std::cin >> y;
//    std::cout  << "Enter z: ";
//    std::cin >> z;
//    std::cout << "Enter c1: ";
//    std::cin >> c1;
//    std::cout << "Enter c2: ";
//    std::cin >> c2;
//    std::cout  << "Enter c3: ";
//    std::cin >> c3;

    
    
    double var[] = {x,y,z};
    std::vector<Operand*> testOperands;
    
    std::vector<Operand*> store;
    std::vector<Operand*> terms;
    
    LeafOperand* c1Coeff = new LeafOperand(c1);
    LeafOperand* c2Coeff = new LeafOperand(c2);
    LeafOperand* c3Coeff = new LeafOperand(c3);

    LeafOperand* x0Var = new LeafOperand(0);
    LeafOperand* x1Var = new LeafOperand(1);
    LeafOperand* x2Var = new LeafOperand(2);
    
    store.push_back(x0Var);
    store.push_back(x1Var);
    AdditionOp* x0p1 = new AdditionOp(2,store);
    store.erase(store.begin(),store.end());
    ExpOp* x0p1e2 = new ExpOp(x0p1,2);
    ExpOp* x0p1e3 = new ExpOp(x0p1,3);
    ExpOp* x0e2 = new ExpOp(x0Var,2);
    
    store.push_back(x0e2);
    store.push_back(x2Var);
    AdditionOp* x0e2p3 = new AdditionOp(2,store);
//    ExpOp* x0e2p3e3 = new ExpOp(x0e2p3,3);
    store.erase(store.begin(),store.end());
    
    
    
    store.push_back(c1Coeff);
    store.push_back(x0p1e2);
    terms.push_back(new MultOp(store.size(),store));
    store.erase(store.begin(),store.end());
    
    store.push_back(c2Coeff);
    store.push_back(x1Var);
    store.push_back(x2Var);
    terms.push_back(new MultOp(store.size(),store));
    store.erase(store.begin(),store.end());
    
    store.push_back(c3Coeff);
    store.push_back(x1Var);
    store.push_back(x0p1e3);
    terms.push_back(new MultOp(store.size(),store));
    store.erase(store.begin(),store.end());
    
    AdditionOp tree(3, terms);
    std::cout << tree.print().str() << std::endl;

    std::cout << "result =  " << tree.evaluate(var) << "\n\n";
    
    
    
    
    Operand* factorSymb = x0p1;
    std::cout << "Factoring out " << factorSymb->print().str() << std::endl;
    std::vector<TreeIndex> index;
    index = Factor::findSymbols(tree, factorSymb);
    
    std::vector<Operand*> treeOperands = tree.getOperands();
    AdditionOp factoredTree;
    
    
    for (int ii = 0; ii < index.size(); ii++)
    {
        Operand* term = treeOperands[index[ii].index];
        store.push_back(Factor::reduceTerm(term, factorSymb));
    }
    
    //Add the terms that were not factored to the factored tree
    for(int ii = 0; ii < treeOperands.size(); ii++)
    {
        bool notFactor = true;
        for(int jj=0; jj < index.size(); jj++)
        {
            if(ii == index[jj].index)
            {
                notFactor = false;
            }
        }
        
        if(notFactor)
        {
            factoredTree.addOperand(treeOperands[ii]);
        }
    }
    
    AdditionOp factoredTerms(store.size(), store);
    store.erase(store.begin(),store.end());
    store.push_back(factorSymb);
    store.push_back(&factoredTerms);
    MultOp* multTerms = new MultOp(2,store);
    store.erase(store.begin(),store.end());
    
    factoredTree.addOperand(multTerms);
    
    std::cout << factoredTree.print().str() << std::endl;
    
    
    
    
    
    
    
 
    
    
    
//    ExpOp* x5Var = new ExpOp(x0Var,5);
//    store.push_back(c1Coeff);
//    store.push_back(x5Var);
//    terms[0] = new MultOp(2,store);
//    testOperands.push_back(terms[0]);
//
//    store.erase(store.begin(),store.end());
//    LeafOperand* yVar = new LeafOperand(1);
//    store.push_back(c2Coeff);
//    store.push_back(yVar);
//    terms[1] = new MultOp(2,store);
//    testOperands.push_back(terms[1]);
//    
//    store.erase(store.begin(),store.end());
//    LeafOperand* zVar = new LeafOperand(2);
//    ExpOp* z2Var = new ExpOp(zVar,2);
//    store.push_back(c3Coeff);
//    store.push_back(z2Var);
//    terms[2] = new MultOp(2,store);
//    testOperands.push_back(terms[2]);
//
//    AdditionOp add(3, testOperands);
//    std::cout << "result =  " << add.evaluate(var) << "\n";
//    
//    std::cout << add.print().str() << std::endl;
    
    
    return 1;
}