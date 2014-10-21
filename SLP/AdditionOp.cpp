//
//  AdditionOp.cpp
//  Bertini2
//
//  Created by James Collins on 11/19/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "LeafOperand.h"
#include "AdditionOp.h"
#include <string>
#include <cstring>
#include <sstream>


/***************  Constructor  ********************
 Input: Vector of Operand pointers
 
 Desc: Stores the input vector as the branches of the Addition node.
 All coefficients are set to 1.
**********************************************/



AdditionOp::AdditionOp(std::vector<Operand*> inOperands)
{
    numOperands = inOperands.size();
    operands = inOperands;
    type = TYPE_ADD;
    //fill coeffs with 1
    for(int ii = 0; ii < numOperands; ii++)
    {
        coeffs.push_back(1);
    }
}





/***************  Constructor  ********************
 Input: Vector of Operand pointers
        Vector of integers, either 1 or -1, representing addition or subtraction
 
 Desc: Stores the input vector as the branches of the Addition node.
 Addition or subtraction is used based on the input vector of integers.
 **********************************************/

AdditionOp::AdditionOp(std::vector<Operand*> inOperands, std::vector<int> inCoeffs)
{
    numOperands = inOperands.size();
    operands = inOperands;
    type = TYPE_ADD;
    coeffs = inCoeffs;
}






/***************  evaluate(Symbol* vars[])  ********************
 Input: Array to Symbol pointers representing values of variables
 
 Output: Symbol pointer to the value of the polynomial
 
 Desc: Evaluates the Addition node recursively.  Calls evaluate on each branch, then 
 adds the results together.
 **********************************************/

Symbol* AdditionOp::evaluate(Symbol* vars[])
{
    Symbol* retval;
    if(operands.size() > 0)
    {
        retval = (operands[0]->evaluate(vars))->clone();
        if(coeffs[0] == -1)
        {
            retval = retval->neg();
        }
        else if(coeffs[0] != 1)
        {
            Symbol* thisCoeff = new IntSymb(coeffs[0]);
            retval = retval->mult(thisCoeff);
        }
    }
    else
    {
        retval = new DoubleSymb(0.0);
    }

    
    for(int ii = 1; ii < operands.size(); ii++)
    {
        if(coeffs[ii] == -1)
        {
            retval = retval->sub(operands[ii]->evaluate(vars));
        }
        else if(coeffs[ii] != 1)
        {
            retval = retval->add(operands[ii]->evaluate(vars)->mult(new IntSymb(coeffs[ii])));
        }
        else
        {
            retval = retval->add(operands[ii]->evaluate(vars));
        }
    }
//    std::cout << std::endl;
    
    return retval;
}






/***************  print()  ********************
 Input: None
 
 Output: stringstream containing character output of polynomial
 
 Desc: Prints the Addition node recursively.  Calls print on each branch, then
 concatenates the strings together with a "+" sign or "-" sign.
 **********************************************/

std::stringstream AdditionOp::print()
{
    std::stringstream ret;
    if(numOperands > 0)
    {
        std::stringstream term = operands[0]->print();
        
        if(coeffs[0] == -1)
        {
            ret << "-" << term.str();
        }
        else if(coeffs[0] != 1)
        {
            ret << coeffs[0] << "("<<term.str() << ")";
        }
        else
        {
            ret << term.str();
        }
        
        
        for(int ii = 1; ii < operands.size(); ii++)
        {
            term = operands[ii]->print();
            if(term.str() != "")
            {
                if(coeffs[ii] == -1)
                {
                    ret << " - " << term.str();
                }
                else if(coeffs[ii] < -1)
                {
                    ret << " - " << abs(coeffs[ii]) << "(" << term.str() << ")";
                }
                else if(coeffs[ii] > 1)
                {
                    ret << " + " << coeffs[ii] << "(" << term.str() << ")";
                }
                else
                {
                    ret << " + " << term.str();
                }
            }
        }
    }
    else
    {
        ret << "";
    }
    
    return ret;
}






/***************  diff(int varIndex)  ********************
 Input: the index of the variable to be differentiated with respect to
 
 Output: Addition node with differentiated polynomial
 
 Desc: Differentiates the Addition node recursively.  Calls diff on each branch, then
 stores each result as a branch of the return Addition node.
 **********************************************/

AdditionOp* AdditionOp::diff(int varIndex)
{
    AdditionOp* retAdd = new AdditionOp();
    
    for(int ii = 0; ii < numOperands; ii++)
    {
        retAdd->addOperand(operands[ii]->diff(varIndex), coeffs[ii]);
    }
    
    
    
    return retAdd;
    
}





/***************  combineAdds()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Addition node.  If the Addition node points to another
 Addition node, method combines branches so that one Addition nodes points to 
 all branches.  Also recursively calls each branch to combine adds for lower down 
 Addition nodes.
 **********************************************/

void AdditionOp::combineAdds()
{
    //Combine Adds for each branch
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->combineAdds();
    }
    
    //Now do the actual work
    int isAdd[numOperands];
    int anyAdds = 0;
    for(int ii = 0; ii < numOperands; ii++)
    {
        if(operands[ii]->getType() == TYPE_ADD)
        {
            isAdd[ii] = 1;
            anyAdds = 1;
        }
        else
        {
            isAdd[ii] = 0;
        }
    }
    
    //Did we find any adds in the branches?
    if(anyAdds)
    {
        //The new branches for the Addition node
        std::vector<Operand*> terms;
        std::vector<int> termCoeffs;
        for(int ii = 0; ii < numOperands; ii++)
        {
            if(isAdd[ii])
            {
                //Add branches for inner Addition node
                AdditionOp* add = (AdditionOp*)operands[ii];
                std::vector<Operand*> innerTerms = operands[ii]->getOperands();
                std::vector<int> innerCoeffs = add->getCoeffs();
                for(int jj = 0; jj < operands[ii]->getNumOperands(); jj++)
                {
                    terms.push_back(innerTerms[jj]);
                    termCoeffs.push_back(innerCoeffs[jj]);
                }
            }
            else
            {
                terms.push_back(operands[ii]);
                termCoeffs.push_back(coeffs[ii]);
            }
        }
        
        //Clear current vector and change to terms
        clear();
        for(int ii = 0; ii < terms.size(); ii++)
        {
            addOperand(terms[ii],termCoeffs[ii]);
        }
    }
    else
    {
        return;
    }
}




/***************  combineMults()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Addition node.  Recursively calls each branch to combine
 mults for nodes lower down the tree.
 Addition nodes. **********************************************/

void AdditionOp::combineMults()
{
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->combineMults();
    }
}








/***************  cleanMults()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Addition node.  Recursively calls each branch to clean
 mults for nodes lower down the tree.
****************************************************/

void AdditionOp::cleanMults()
{
    int findZero = 0;
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->cleanMults();
        if(operands[ii]->getType() == TYPE_MULT)
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
        cleanAdds();
    }
}











/***************  cleanAdds()  ********************
 Input: None
 
 Output: None
 
 Desc: Modifies the Addition node.  Recursively calls each branch to clean
 additions for nodes lower down the tree.  If there is a 0 term, entire addition node
 points to a 0 leaf.
 ****************************************************/

void AdditionOp::cleanAdds()
{
    //clean Mults for each branch.  If we clean a branch to 0, erase that branch
    for(int ii = 0; ii < numOperands; ii++)
    {
        operands[ii]->cleanAdds();
    }
    
    
    
    
    
    //If there is a 0 as a term, erase it
    std::vector<Operand*> store;
    std::vector<int> storeCoeff;
    for(int ii = 0; ii < numOperands; ii++)
    {
        int isZero = 0;
        if(operands[ii]->getType() == TYPE_LEAF)
        {
            LeafOperand* leaf = (LeafOperand*)operands[ii];
            if(leaf->getIsZero())
            {
                isZero = 1;
            }
        }
        //If it isn't a 0, store in vector
        if(!isZero)
        {
            store.push_back(operands[ii]);
            storeCoeff.push_back(coeffs[ii]);
        }
    }
    
    //clear operands and replace with store
    clear();
    if(store.size() > 0)
    {
        for(int ii = 0; ii < store.size(); ii++)
        {
            addOperand(store[ii],storeCoeff[ii]);
        }
    }
    else
    {
        addOperand(new LeafOperand(0,1),1);
    }
    
}






