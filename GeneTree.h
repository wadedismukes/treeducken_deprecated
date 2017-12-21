//
//  GeneTree.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef GeneTree_h
#define GeneTree_h

#include "LocusTree.h"

class GeneTree : public Tree {
    private:
        unsigned numTips;
        unsigned individualsPerPop;
        unsigned popSize;
    
    public:
        GeneTree(MbRandom *rando, unsigned nt, unsigned ipp, unsigned ne);
        ~GeneTree();
};

#endif /* GeneTree_hpp */
