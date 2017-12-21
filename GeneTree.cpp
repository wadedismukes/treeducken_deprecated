//
//  GeneTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "GeneTree.h"

GeneTree::GeneTree(MbRandom *p, unsigned nt, unsigned ipp, unsigned ne) : Tree(p, nt){
    rando = p;
    numTaxa = nt;
    individualsPerPop = ipp;
    popSize = ne;
    
}

GeneTree::~GeneTree(){
    
}

