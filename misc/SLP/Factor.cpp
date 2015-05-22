//
//  Factor.cpp
//  Bertini2
//
//  Created by James Collins on 12/10/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#include "Factor.h"


//Only for test
Operand* Factor::testFactor = 0;











//Determine whether a factor is within a particular term
bool Factor::isFactorInTerm(Operand* term, Operand* factor)
{
    bool isFactorTerm = false;
    
    Operand* thisTerm = term;
    
    //Is symb in thisTerm?
    if(thisTerm->getType() == Operand::TYPE_LEAF)
    {
        if(thisTerm == factor)
        {
            isFactorTerm = true;
        }
    }
    else if(thisTerm->getType() == Operand::TYPE_EXP)
    {
        ExpOp* thisExpTerm = (ExpOp*)thisTerm;
        if(thisExpTerm->getBase() == factor)
        {
            isFactorTerm = true;
        }
    }
    else
    {
        //loop thru factors of the mult term
        std::vector<Operand*> factors = thisTerm->getOperands();
        for(int jj = 0; jj < factors.size(); jj++)
        {
            Operand* thisFactor = factors[jj];
            if(thisFactor->getType() == Operand::TYPE_LEAF)
            {
                if(thisFactor == factor)
                {
                    isFactorTerm = true;
                }
            }
            else if(thisFactor->getType() == Operand::TYPE_EXP)
            {
                ExpOp* thisExpFactor = (ExpOp*)thisFactor;
                if(thisExpFactor->getBase() == factor)
                {
                    isFactorTerm = true;
                }
            }
        }
    }
    
    
    
    return isFactorTerm;
}





//Determine whether a factor is within a particular term, if so how many times?
bool Factor::isFactorInTerm(Operand* term, Operand* factor, int &expCount)
{
    bool isFactorTerm = false;
    expCount = -1;
    
    Operand* thisTerm = term;
    
    //Is symb in thisTerm?
    if(thisTerm->getType() == Operand::TYPE_LEAF)
    {
        if(thisTerm == factor)
        {
            isFactorTerm = true;
            expCount = 1;
        }
    }
    else if(thisTerm->getType() == Operand::TYPE_EXP)
    {
        ExpOp* thisExpTerm = (ExpOp*)thisTerm;
        if(thisExpTerm->getBase() == factor)
        {
            isFactorTerm = true;
            expCount = thisExpTerm->getExp();
        }
    }
    else
    {
        //loop thru factors of the mult term
        std::vector<Operand*> factors = thisTerm->getOperands();
        for(int jj = 0; jj < factors.size(); jj++)
        {
            Operand* thisFactor = factors[jj];
            if(thisFactor->getType() == Operand::TYPE_LEAF)
            {
                if(thisFactor == factor)
                {
                    isFactorTerm = true;
                    expCount = 1;
                }
            }
            else if(thisFactor->getType() == Operand::TYPE_EXP)
            {
                ExpOp* thisExpFactor = (ExpOp*)thisFactor;
                if(thisExpFactor->getBase() == factor)
                {
                    isFactorTerm = true;
                    expCount = thisExpFactor->getExp();
                }
            }
        }
    }
    
    
    
    return isFactorTerm;
}








//Find all symbols in an AdditionOp
//Output is an int array whose first index is the symbol and second is the term
std::vector<Operand*> Factor::findAllSymbols(AdditionOp* node)
{
    /* Initialize variables */
    //N = number of terms
    int N = node->getNumOperands();
    //Symbols = list of distinct symbols in all terms
    std::vector<Operand*> symbols;
    
    
    
    std::vector<Operand*> terms = node->getOperands();
    for(int ii = 0; ii < N; ii++)
    {
        //If term is just a leaf(single factor)
        if(terms[ii]->getType() == Operand::TYPE_LEAF)
        {
            //Check if this leaf is in symbols
            bool inSymbols = false;
            for(int jj = 0; jj < symbols.size(); jj++)
            {
                if(terms[ii] == symbols[jj])
                {
                    inSymbols = true;
                }
            }
            if(!inSymbols)
            {
                symbols.push_back(terms[ii]);
            }
        }
        
        //If term is an ExpOp(single factor raised to an exponent)
        else if(terms[ii]->getType() == Operand::TYPE_EXP)
        {
            //Check if the base is in symbols
            bool inSymbols = false;
            ExpOp* expTerm = (ExpOp*)terms[ii];
            for(int jj = 0; jj < symbols.size(); jj++)
            {
                if(expTerm->getBase() == symbols[jj])
                {
                    inSymbols = true;
                }
            }
            if(!inSymbols)
            {
                symbols.push_back(expTerm->getBase());
            }
        }
        
        
        //If term is an MultOp(many factors)
        else if(terms[ii]->getType() == Operand::TYPE_MULT)
        {
            
            std::vector<Operand*> factors = terms[ii]->getOperands();
            for(int kk = 0; kk < factors.size(); kk++)
            {
                //Check if the factor is in symbols
                bool inSymbols = false;
                Operand* thisFactor;
                if(factors[kk]->getType() == Operand::TYPE_EXP)
                {
                    ExpOp* expFactor = (ExpOp*)factors[kk];
                    thisFactor = expFactor->getBase();
                }
                else
                {
                    thisFactor = factors[kk];
                }
                for(int jj = 0; jj < symbols.size(); jj++)
                {
                    if(thisFactor == symbols[jj])
                    {
                        inSymbols = true;
                    }
                    

                } //end symbols loop
                if(!inSymbols)
                {
                    symbols.push_back(thisFactor);
                }

            } //end factors loop
        } //endif mult type
    } //end terms loop
    
    return symbols;
}








//Count how many times each symbols appear in each term and return array
std::array<matrix<int>,2> Factor::countSymbols(AdditionOp* node, std::vector<Operand*> symbols)
{
    int N = node->getNumOperands();
    int M = symbols.size();
    int expCount;
    
    if(N == 0 || M == 0)
    {
        std::array<matrix<int>,2> ret;
        return ret;
    }
    
    std::array<matrix<int>,2> retArray = {matrix<int>(N,M), matrix<int>(N,M)};
    std::vector<Operand*> terms = node->getOperands();
    
    
    //Loop thru all the terms
    for(int ii = 0; ii < N; ii++)
    {
        //Loop thru all the symbols
        for(int jj = 0; jj < M; jj++)
        {
            //If the symbol is in the term...
            if(isFactorInTerm(terms[ii], symbols[jj], expCount))
            {
                //...store it's exponent in array
                retArray[0](ii,jj) = expCount;
                retArray[1](ii,jj) = 1;
            }
            else
            {
                retArray[0](ii,jj) = 0;
                retArray[1](ii,jj) = 0;
            }
        }//Endl symbols loop
    }//End terms loop
    
    retArray[0] = trans(retArray[0]);
    retArray[1] = trans(retArray[1]);
    
    return retArray;
}







AdditionOp* Factor::factorAdd(AdditionOp* node)
{
    //Determine what to factor
    Operand* factor = findFactor(node);
    
    //If not factor, return node
    if(factor == 0)
    {
        return node;
    }
    
    // Split the terms into terms with factor and terms without
    AdditionOp* factoredTerms = new AdditionOp();
    AdditionOp* restTerms = new AdditionOp();
    splitTree(node, factor, factoredTerms, restTerms);
    
    //Now recursively apply this algorithm to factoredTerms terms
    factoredTerms = factorAdd(factoredTerms);
    
    //Create MultOp
    MultOp* factorMult = new MultOp();
    factorMult->addOperand(factor);
    factorMult->addOperand(factoredTerms);
    
    //Finally we deal with restTerms...
    
    //...if there are more than 1 term left...
    if(restTerms->getNumOperands() > 1)
    {
        //Recursively apply algorithm to factor the rest of the terms
        restTerms = factorAdd(restTerms);
    }
    
    AdditionOp* newNode = new AdditionOp();
    newNode->addOperand(factorMult);
    newNode->addOperand(restTerms);

    
    
    
    return newNode;
}







//Determine what to factor out of a set of terms
Operand* Factor::findFactor(AdditionOp* baseAdd)
{
    std::vector<Operand*> symbols = findAllSymbols(baseAdd);
    std::array<matrix<int>,2> countMatrix = countSymbols(baseAdd, symbols);
    
    int N = countMatrix[0].size1();
    int M = countMatrix[0].size2();
    vector<int> ones(M);
    for(int ii = 0; ii < ones.size(); ii++)
    {
        ones(ii) = 1;
    }
    
    //Calculate symbol sums...
    vector<int> sumExp = prod(countMatrix[0],ones);
    vector<int> sumTerms = prod(countMatrix[1],ones);
    //Is there a factor?
    int max = 0;
    for(int ii = 0; ii < sumTerms.size(); ii++)
    {
        if(sumTerms(ii) > max)
        {
            max = sumTerms(ii);
        }
    }
    if(max < 2)
    {
        //There is nothing to factor out of terms
        return 0;
    }
    vector<int> sumSymbols = sumExp + sumTerms;
    
    sumSymbols = sumTerms;
    
    //...find which symbol has the max...
    max = 0; int maxInd = -1;
    for(int ii = 0; ii < sumSymbols.size(); ii++)
    {
        if(sumSymbols(ii) > max)
        {
            max = sumSymbols(ii);
            maxInd = ii;
        }
    }
    
    return symbols[maxInd];
    
    
    
    
    
    
    
    
//    int counter = 0;
//    //Only for test
//    for(int ii = 0; ii < baseAdd->getNumOperands(); ii++)
//    {
//        if(isFactorInTerm(baseAdd->getOperands()[ii], Factor::testFactor))
//        {
//            counter++;
//        }
//    }
//    
//    if(counter >= 2)
//    {
//        return Factor::testFactor;
//    }
//    else
//    {
//        return 0;
//    }
    
    
}




//Split node leaves into factored terms(factorAdd) and the rest(restAdd)
void Factor::splitTree(AdditionOp* baseAdd, Operand* symb, AdditionOp* factorAdd, AdditionOp* restAdd)
{
    //loop thru all terms in baseAdd
    std::vector<Operand*> terms = baseAdd->getOperands();
    std::vector<int> coeffs = baseAdd->getCoeffs();
    for(int ii = 0; ii < terms.size(); ii++)
    {
        Operand* thisTerm = terms[ii];
        
        
        //Is symb in thisTerm?
        bool isFactorTerm = isFactorInTerm(thisTerm, symb);
        
        //If term has factor in it...
        if(isFactorTerm)
        {
            
            //...divide the term by the factor..
            Operand* factoredTerm = reduceTerm(thisTerm, symb);
            //..and push the now factored term onto factorAdd
            factorAdd->addOperand(factoredTerm, coeffs[ii]);
        }
        //If term does not have factor in it...
        else
        {
            //...push it onto restAdd
            restAdd->addOperand(thisTerm, coeffs[ii]);
        }
        
    }
    
}







//Only works if symb is a single factor (i.e. x, not x^2)
Operand* Factor::reduceTerm(Operand* term, Operand* symb)
{
    MultOp* ret = new MultOp();
    
    
    //If leaf no factoring to be done
    if(term->getType() == Operand::TYPE_LEAF)
    {
        if(term == symb)
        {
            return new LeafOperand(new IntSymb(1));
        }
        else
        {
            return term;
        }
    }
    //If exp...
    else if(term->getType() == Operand::TYPE_EXP)
    {
        ExpOp* expTerm = (ExpOp*)term;
        if(isFactorInTerm(expTerm, symb))
        {
            int expNow = expTerm->getExp();
            if(expNow > 2)
            {
                ExpOp* newExp = new ExpOp(expTerm->getBase(), expNow-1);
                return newExp;
            }
            else if(expNow == 2)
            {
                Operand* leafTerm = expTerm->getBase();
                return leafTerm;
            }
            else
            {
                return new LeafOperand(new IntSymb(1));
            }
        }
    }
    //Otherwise...
    else
    {
        std::vector<Operand*> ops = term->getOperands();
        for(int ii = 0; ii < ops.size(); ii++)
        {
            if(ops[ii]->getType() != Operand::TYPE_EXP)
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
                        ExpOp* newExp = new ExpOp(exop->getBase(), exp-1);
                        ret->addOperand(newExp);
                    }
                }
                else
                {
                    ret->addOperand(exop);
                }
            }
        }
        
        //If reduced term is 1
        if(ret->getNumOperands() == 0)
        {
            return new LeafOperand(new IntSymb(1));
        }
    }
    
    
    
    
    
    
    return ret;
}




int Factor::countAdds(Operand* node)
{
    int type = node->getType();
    
    if(type == Operand::TYPE_LEAF)
    {
        return 0;
    }
    else if(type == Operand::TYPE_EXP)
    {
        ExpOp* expNode = (ExpOp*)node;
        return countAdds(expNode->getBase());
    }
    else if(type == Operand::TYPE_MULT)
    {
        int addCount = 0;
        std::vector<Operand*> nodeOps = node->getOperands();
        for(int ii = 0; ii < node->getNumOperands(); ii++)
        {
            addCount += countAdds(nodeOps[ii]);
        }
        return addCount;
    }
    else
    {
        int addCount = -1;
        std::vector<Operand*> nodeOps = node->getOperands();
        for(int ii = 0; ii < node->getNumOperands(); ii++)
        {
            addCount += countAdds(nodeOps[ii]);
            if(nodeOps[ii]->getOperands().size() > 0 || nodeOps[ii]->getType() == Operand::TYPE_LEAF
               || nodeOps[ii]->getType() == Operand::TYPE_EXP)
            {
                addCount += 1;
            }
        }
        return std::max(addCount, 0);
    }
}







int Factor::countMults(Operand* node)
{
    int type = node->getType();
    
    if(type == Operand::TYPE_LEAF)
    {
        return 0;
    }
    else if(type == Operand::TYPE_EXP)
    {
        int multCount = 0;
        ExpOp* expNode = (ExpOp*)node;
        multCount += countAdds(expNode->getBase());
        if(expNode->getBase() != 0)
        {
            multCount += expNode->getExp() - 1;
        }
        
        return multCount;
    }
    else if(type == Operand::TYPE_MULT)
    {
        int multCount = -1;
        //Get the multiplication operands
        std::vector<Operand*> nodeOps = node->getOperands();
        for(int ii = 0; ii < node->getNumOperands(); ii++)
        {
            multCount += countMults(nodeOps[ii]);
            if(nodeOps[ii]->getOperands().size() > 0 || nodeOps[ii]->getType() == Operand::TYPE_LEAF
               || nodeOps[ii]->getType() == Operand::TYPE_EXP)
            {
                multCount += 1;
            }
        }
        return std::max(multCount, 0);
    }
    else
    {
        int multCount = 0;
        std::vector<Operand*> nodeOps = node->getOperands();
        for(int ii = 0; ii < node->getNumOperands(); ii++)
        {
            multCount += countMults(nodeOps[ii]);
        }
        return multCount;

    }
}

