//
//  TreeIndex.h
//  Bertini2
//
//  Created by James Collins on 12/10/13.
//  Copyright (c) 2013 James Collins. All rights reserved.
//

#ifndef Bertini2_TreeIndex_h
#define Bertini2_TreeIndex_h

class TreeIndex
{
public:
    TreeIndex() {};
    TreeIndex(int inIndex, bool inIsExp) {index = inIndex; isExp = inIsExp;};
    int index; //Represents the index of the term in the tree where index is found
    bool isExp; //Is the symbol found in an exponent operator or alone?
    int mult; // Multiplicity of the symbol at that index
    
};


#endif
