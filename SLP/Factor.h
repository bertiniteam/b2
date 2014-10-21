//
//  Factor.h
//  Bertini2
//
//  Created by James Collins on 12/10/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef __Bertini2__Factor__
#define __Bertini2__Factor__

#include <iostream>
#include <vector>
#include <array>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Operand.h"
#include "AdditionOp.h"
#include "MultOp.h"
#include "ExpOp.h"
#include "LeafOperand.h"
#include "Symbol.h"
#include "IntSymb.h"

using namespace boost::numeric::ublas;


//All of Factor assumes that the base node of the poly tree is Add and the immediate leaves of the base node are either Mult, Exp or a Leaf.
//If this is not the case, then the poly tree can easily be put in this form

class Factor
{
private:
    //Determine if factor is in a particular term
    static bool isFactorInTerm(Operand* term, Operand* factor);
    
    //Determine if factor is in a particular term, if so how many times?
    static bool isFactorInTerm(Operand* term, Operand* factor, int &expCount);
    
    //Find all symbols in the terms in node
    static std::vector<Operand*> findAllSymbols(AdditionOp* node);
    
    //Count how many times each symbols appear in each term and return array
    static std::array<matrix<int>,2> countSymbols(AdditionOp* node, std::vector<Operand*> symbols);
    
public:
    
    
    
    //Remove factor from term
    static Operand* reduceTerm(Operand* term, Operand* symb);
    
    
    //The recursive algorithm that factors an addition op
    static AdditionOp* factorAdd(AdditionOp* node);

    
    //Determine what to factor out of a set of terms
    static Operand* findFactor(AdditionOp* terms);
    
    //Split the tree into a part that can factor out symb and a part that cannot
    static void splitTree(AdditionOp* baseAdd, Operand* symb, AdditionOp* factorAdd, AdditionOp* restAdd);
    
    //Count the additions in a tree
    static int countAdds(Operand* node);
    
    //Count the multiplications in a tree
    static int countMults(Operand* node);
    
    
    //Only for the test
    static Operand* testFactor;
};



#endif /* defined(__Bertini2__Factor__) */
