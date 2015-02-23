//
//  Factor.cpp
//  Bertini2
//
//  Created by James Collins on 12/10/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "Factor.h"


std::vector<TreeIndex> Factor::findSymbols(AdditionOp tree, Operand* symb)
{
    std::vector<TreeIndex> ret;
    std::vector<Operand*> terms = tree.getOperands();
        
    for(int ii = 0; ii < tree.getNumOperands(); ii++)
    {
        std::vector<Operand*> operands = terms[ii]->getOperands();
        for(int jj = 0; jj < terms[ii]->getNumOperands(); jj++)
        {
            if(operands[jj]->IsExp())
            {
                ExpOp* op = (ExpOp*)operands[jj];
                Operand* base = op->getBase();
                if(base == symb)
                {
                    TreeIndex index(ii,true);
                    ret.push_back(index);
                    break;
                }
            }
            else
            {
                if(operands[jj] == symb)
                {
                    TreeIndex index(ii,false);
                    ret.push_back(index);
                    break;
                }
            }
        }
    }
    
    return ret;
}





Operand* Factor::reduceTerm(Operand* term, Operand* symb)
{
    MultOp* ret = new MultOp();
    
    
    if(term->IsLeaf())
    {
        ret = 0;
        return ret;
    }
    else
    {
        std::vector<Operand*> ops = term->getOperands();
        for(int ii = 0; ii < ops.size(); ii++)
        {
            if(!ops[ii]->IsExp())
            {
                if(ops[ii] != symb)
                {
                    ret->addOperand(ops[ii]);
                }
            }
            else
            {
                ExpOp* exop = (ExpOp*)ops[ii];
                if(exop->getBase() == symb)
                {
                    int exp = exop->getExp();
                    if(exp-1 == 1)
                    {
                        ret->addOperand(exop->getBase());
                    }
                    else
                    {
                        exop->setExp(exp-1);
                        ret->addOperand(exop);
                    }
                }
                else
                {
                    ret->addOperand(exop);
                }
            }
        }
    }
    
    
    //if we just have a leaf
    std::vector<Operand*> ops = ret->getOperands();
    if(ops.size() == 1)
    {
        delete ret;
        Operand* ret = ops[0];
        return ret;
    }
    
    if(ops.size() == 0)
    {
        delete ret;
        Operand* ret = 0;
        return ret;
    }
    
    
    
    
    return ret;
}





void Factor::splitTree(AdditionOp baseTree, Operand* symb, AdditionOp* factorTree, AdditionOp* restTree)
{
    
}















