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
#include "Operand.h"
#include "AdditionOp.h"
#include "MultOp.h"
#include "ExpOp.h"
#include "LeafOperand.h"
#include "TreeIndex.h"

class Factor
{
public:
    //Go thru tree and find all terms with the symbol
    static std::vector<TreeIndex> findSymbols(AdditionOp tree, Operand* symb);
    
    //Remove factor from term
    static Operand* reduceTerm(Operand* term, Operand* symb);
    
    //Split the tree into a part that can factor out symb and a part that cannot
    static void splitTree(AdditionOp baseTree, Operand* symb, AdditionOp* factorTree, AdditionOp* restTree);
    
    //
};

#endif /* defined(__Bertini2__Factor__) */
