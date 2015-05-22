//
//  MultOp.cpp
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "MultOp.h"
#include "LeafOperand.h"



/***************  Constructor  ********************
 Input: Vector of Operand pointers
 
 Desc: Stores the input vector as the branches of the Mult node.
 **********************************************/

MultOp::MultOp(std::vector<Operand*> inOperands)
{
    numOperands = inOperands.size();
    operands = inOperands;
    type = TYPE_MULT;
}










/***************  evaluate(Symbol* vars[])  ********************
 Input: Array to Symbol pointers representing values of variables
 
 Output: Symbol pointer to the value of the polynomial
 
 Desc: Evaluates the Mult node recursively.  Calls evaluate on each branch, then
 multiplies the results together.
 **********************************************/

Symbol* MultOp::evaluate(Symbol* vars[])
{
    Symbol* retval;
    if(operands.size() > 0)
    {
        retval = operands[0]->evaluate(vars)->clone();
    }
    else{
        retval = new DoubleSymb(1.0);
    }
//    std::cout << "Mult test: " << operands[0]->evaluate(vars)->print().str() << " ";
    
    for(int ii = 1; ii < operands.size(); ii++)
    {
        retval = retval->mult(operands[ii]->evaluate(vars));
//        std::cout << operands[ii]->evaluate(vars)->print().str() << " ";
    }
    
    return retval;
}









/***************  print()  ********************
 Input: None
 
 Output: stringstream containing character output of polynomial
 
 Desc: Prints the Mult node recursively.  Calls print on each branch, then
 concatenates the strings together with a "*" sign.
 **********************************************/

std::stringstream MultOp::print()
{
    std::stringstream ret;
    if(operands.size() > 0)
    {
        std::stringstream term = operands[0]->print();
        
        if(operands.size() == 1)
        {
            ret << term.str();
        }
        else
        {
            if(operands[0]->getType() == TYPE_LEAF)
            {
                ret << term.str();
            }
            else
            {
                ret << "(" << term.str() << ")";
            }
            
            for(int ii = 1; ii < operands.size(); ii++)
            {
                term = operands[ii]->print();
                if(operands[ii]->getType() == TYPE_LEAF)
                {
                    ret << term.str();
                }
                else
                {
                    ret << "(" << term.str() << ")";
                }
            }
 
        }
    }
    
    return ret;
}










/***************  diff(int varIndex)  ********************
 Input: the index of the variable to be differentiated with respect to
 
 Output: Addition node with differentiated polynomial
 
 Desc: Differentiates the Mult node recursively.  Calls diff on each branch, then
 combines them using the Leibnitz rule.
 **********************************************/

AdditionOp* MultOp::diff(int varIndex)
{
    AdditionOp* retAdd = new AdditionOp();
    
    for(int ii = 0; ii < numOperands; ii++)
    {
        MultOp* multTemp = new MultOp();
        Operand* diff = operands[ii]->diff(varIndex);
        if(diff->getType() == TYPE_LEAF)
        {
            LeafOperand* leaf = (LeafOperand*)diff;
            if(leaf->getIsZero())
            {
                continue;
            }
        }
        for(int jj = 0; jj < numOperands; jj++)
        {
            if(jj != ii)
            {
                multTemp->addOperand(operands[jj]);
            }
        }
        multTemp->addOperand(diff);
        
        retAdd->addOperand(multTemp);
    }
    
    
    
    
    return retAdd;
    
}










/***************  combineAdds()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Mult node.  Recursively calls each branch to combine
 adds for nodes lower down the tree.
 **********************************************/

void MultOp::combineAdds()
{
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->combineAdds();
    }
}













/***************  combineMults()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Mult node.  If the Mult node points to another
 Mult node, method combines branches so that one Mult nodes points to
 all branches.  Also recursively calls each branch to combine mults for lower down
 Mult nodes.
******************************************************/

void MultOp::combineMults()
{
    //Combine Mults for each branch
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->combineAdds();
    }
    
    //Now do the actual work
    int isMult[numOperands];
    int anyMults = 0;
    for(int ii = 0; ii < numOperands; ii++)
    {
        if(operands[ii]->getType() == TYPE_MULT)
        {
            isMult[ii] = 1;
            anyMults = 1;
        }
        else
        {
            isMult[ii] = 0;
        }
    }
    
    //Did we find any adds in the branches?
    if(anyMults)
    {
        std::vector<Operand*> factors;
        for(int ii = 0; ii < numOperands; ii++)
        {
            if(isMult[ii])
            {
                std::vector<Operand*> innerFactors = operands[ii]->getOperands();
                for(int jj = 0; jj < operands[ii]->getNumOperands(); jj++)
                {
                    factors.push_back(innerFactors[jj]);
                }
            }
            else
            {
                factors.push_back(operands[ii]);
            }
        }
        
        //Clear current vector and change to terms
        clear();
        for(int ii = 0; ii < factors.size(); ii++)
        {
            addOperand(factors[ii]);
        }
    }
    else
    {
        return;
    }
}







/***************  cleanMults()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Mult node.  Recursively calls each branch to clean
 mults for nodes lower down the tree.  If there is a 0 factor, entire mult node 
 points to a 0 leaf.  If there is a 1 factor, erase that branch, unless it is 
 the only branch.
 ****************************************************/

void MultOp::cleanMults()
{
    //clean Mults for each branch.  If we clean a branch to 0, erase that branch
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->cleanMults();
    }
    
    
    //First check if there are any 0's as a factor, if so replace whole thing with 0
    for(int ii = 0; ii < numOperands; ii++)
    {
        if(operands[ii]->getType() == TYPE_LEAF)
        {
            LeafOperand* leaf = (LeafOperand*)operands[ii];
            if(leaf->getIsZero())
            {
                clear();
                addOperand(new LeafOperand(0,1));
                return;
            }
        }
    }
    
    
    
    //If there is a 1 as a factor, erase it
    std::vector<Operand*> store;
    for(int ii = 0; ii < numOperands; ii++)
    {
        int isOne = 0;
        if(operands[ii]->getType() == TYPE_LEAF)
        {
            LeafOperand* leaf = (LeafOperand*)operands[ii];
            if(leaf->getIsOne())
            {
                isOne = 1;
            }
        }
        //If it isn't a 1, store in vector
        if(!isOne)
        {
            store.push_back(operands[ii]);
        }
    }
    
    //clear operands and replace with store
    clear();
    if(store.size() > 0)
    {
        for(int ii = 0; ii < store.size(); ii++)
        {
            addOperand(store[ii]);
        }
    }
    else
    {
        addOperand(new LeafOperand(1,0));
    }
    
}













/***************  cleanAdds()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Mult node.  Recursively calls each branch to clean
 adds for nodes lower down the tree.
 ****************************************************/

void MultOp::cleanAdds()
{
    int findZero = 0;
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->cleanAdds();
        if(operands[ii]->getType() == TYPE_ADD)
        {
            if(operands[ii]->getNumOperands() == 1)
            {
                if(operands[ii]->getOperands()[0]->getType() == TYPE_LEAF)
                {
                    LeafOperand* leaf = (LeafOperand*)operands[ii]->getOperands()[0];
                    if(leaf->getIsZero())
                    {
                        findZero = 1;
                        delete operands[ii];
                        operands[ii] = new LeafOperand(0,1);
                    }
                }
            }
        }
    }
    
    if(findZero)
    {
        cleanMults();
    }
}






