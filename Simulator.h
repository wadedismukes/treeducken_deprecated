//
//  Simulator.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Simulator_h
#define Simulator_h
#include "SpeciesTree.h"
#include "LocusTree.h"
#include <set>
#include <map>

class Simulator
{
protected:
    MbRandom    *rando;
    double      currentSimTime;
    unsigned    simType;
    unsigned    numTaxaToSim, gsaStop;
    unsigned    numLoci;
    unsigned    numGenes;
    double      speciationRate, extinctionRate;
    double      samplingRate;
    double      geneBirthRate, geneDeathRate, transferRate;
    double      propTransfer, propDuplicate;
    unsigned    indPerPop;
    unsigned    popSize;
    std::vector<SpeciesTree*>   gsaTrees;
    SpeciesTree*    spTree;
    std::vector<LocusTree*>     lociTrees;

public:
    // Simulating species tree only
    Simulator(MbRandom *p, unsigned numTaxaToSim, double speciationRate, double extinctionRate, double rho);
    // Simulating species and locus tree
    Simulator(MbRandom *p, unsigned numTaxaToSim, double speciationRate, double extinctionRate, double rho, unsigned numLociToSim, double geneBirthRate, double geneDeathRate, double transferRate);
    // Simulating species and locus tree with proportion of transfer (e.g. hybridization, linkage)
    Simulator(MbRandom *p, unsigned numTaxaToSim, double speciationRate, double extinctionRate, double rho, unsigned numLociToSim, double geneBirthRate, double geneDeathRate, double transferRate, double propTransfer);

    ~Simulator();
    bool    gsaBDSim();
    bool    bdsaBDSim();
    bool    simSpeciesTree();
    bool    simSpeciesLociTrees();
    bool    gsaCheckStop();
    void    initializeSim();
    void    processGSASim();
    void    prepGSATreeForReconstruction();
    void    processSpTreeSim();
    std::string    printSpeciesTreeNewick();
    std::string    printLocusTreeNewick();
};


#endif /* Simulator_h */
