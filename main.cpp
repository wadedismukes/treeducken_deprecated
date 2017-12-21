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
    
    unsigned numTaxa = 100;
    std::string sptree;
    std::string locusTree;
    double br = 2.0;
    double dr = 1.0;
    double gbr = 2.0;
    double gdr = 0.25;
    double lgtr = 1.0;
//    Simulator *sim = new Simulator(&rand, numTaxa, br, dr, 1.0);
//    sim->simSpeciesTree();
//    sptree = sim->printSpeciesTreeNewick();
//    std::cout << "species tree: " << sptree << std::endl;
    Simulator *sim = new Simulator(&rand, numTaxa, br, dr, 1.0, 1, gbr, gdr, lgtr);
    sim->simSpeciesLociTrees();
    sptree = sim->printSpeciesTreeNewick();
    std::cout << "species tree: " << std::endl << sptree << std::endl <<  std::endl;
    locusTree = sim->printLocusTreeNewick();
    std::cout << "Locus Tree: " << std::endl << locusTree << std::endl << std::endl;
    return 0;
}
