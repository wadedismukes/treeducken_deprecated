//
//  main.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/7/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include <iostream>
#include "SpeciesTree.h"
#include "Simulator.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    MbRandom rand;
    int sd1 = 38013;
    int sd2 = 23118;
    if(sd1 > 0 && sd2 > 0)
        rand.setSeed(sd1, sd2);
    else
        rand.setSeed();
    seedType gs1, gs2;
    rand.getSeed(gs1, gs2);
    std::cout << "\nSeeds {" << gs1 << ", " << gs2 << "}" << std::endl;
    
    unsigned numTaxa = 5;
    std::string sptree;
    std::string locusTree;
    std::string geneTree;
    double br = 4.0;
    double dr = 2.0;
    double gbr = 2.0;
    double gdr = 0.25;
    double lgtr = 1.0;
    unsigned ipp = 4;
    unsigned popSize = 10;
    double genTime = 100;
    
    Simulator *sim = new Simulator(&rand, numTaxa, br, dr, 1.0, 1, gbr, gdr, lgtr, ipp, popSize, genTime);
    sim->simThreeTree();
    sptree = sim->printSpeciesTreeNewick();
    std::cout << "species tree: " << std::endl << sptree << std::endl <<  std::endl;
    locusTree = sim->printLocusTreeNewick();
    std::cout << "Locus Tree: " << std::endl << locusTree << std::endl << std::endl;
    geneTree = sim->printGeneTreeNewick();
        std::cout << "Gene Tree: " << std::endl << geneTree << std::endl << std::endl;
    
    return 0;
}
