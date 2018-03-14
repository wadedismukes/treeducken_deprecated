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
    int sd1 = 0;
    int sd2 = 0;
    if(sd1 > 0 && sd2 > 0)
        rand.setSeed(sd1, sd2);
    else
        rand.setSeed();
    seedType gs1, gs2;
    rand.getSeed(gs1, gs2);
    std::cout << "\nSeeds {" << gs1 << ", " << gs2 << "}" << std::endl;
    
    unsigned numTaxa = 4;
    std::string sptree;
    std::string locusTree;
    std::string geneTree;
    double br = 0.01;
    double dr = 0.005;
    double gbr = 0.02;
    double gdr = 0.01;
    double lgtr = 0.02;
    unsigned ipp = 4;
    unsigned popSize = 10000;
    double genTime = 1;
    unsigned numLoci = 10;
    
    Simulator *sim = new Simulator(&rand, numTaxa, br, dr, 1.0, numLoci, gbr, gdr, lgtr, ipp, popSize, genTime);
    sim->simThreeTree();
    sptree = sim->printSpeciesTreeNewick();
    std::cout << "species tree: " << std::endl << sptree << std::endl <<  std::endl;
    for(int i = 0; i < numLoci; i++){
        locusTree = sim->printLocusTreeNewick(i);
        std::cout << "Locus Tree: " << std::endl << locusTree << std::endl << std::endl;
        geneTree = sim->printGeneTreeNewick(i);
        std::cout << "Gene Tree: " << std::endl << geneTree << std::endl << std::endl;
    }
    return 0;
}
