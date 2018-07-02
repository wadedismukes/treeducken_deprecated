//
//  Simulator.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Simulator_h
#define Simulator_h
#include "GeneTree.h"
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
        double      treeScale;
        double      geneBirthRate, geneDeathRate, transferRate;
        double      propTransfer, propDuplicate;
        unsigned    indPerPop;
        unsigned    popSize;
        double      generationTime;
        double      outgroupFrac;
        std::vector<SpeciesTree*>   gsaTrees;
        SpeciesTree*    spTree;
        LocusTree*      lociTree;
        std::vector<LocusTree*> locusTrees;
        GeneTree*       geneTree;
        std::vector<std::vector<GeneTree*> > geneTrees;

    public:
        // Simulating species tree only
        Simulator(MbRandom *p, 
                  unsigned numTaxaToSim, 
                  double speciationRate, 
                  double extinctionRate, 
                  double rho);
        // Simulating species and locus tree
        Simulator(MbRandom *p, 
                  unsigned numTaxaToSim, 
                  double speciationRate, 
                  double extinctionRate, 
                  double rho, 
                  unsigned numLociToSim, 
                  double geneBirthRate, 
                  double geneDeathRate, 
                  double transferRate);
        // Simulating species and locus tree with proportion of transfer (e.g. hybridization, linkage)
        Simulator(MbRandom *p,
                unsigned numTaxaToSim,
                double speciationRate,
                double extinctionRate,
                double rho, 
                unsigned numLociToSim, 
                double geneBirthRate, 
                double geneDeathRate, 
                double transferRate, 
                double propTransfer);
        // Simulating species and locus trees with one gene tree per locus tree
        Simulator(MbRandom *p,
                unsigned numTaxaToSim,
                double speciationRate,
                double extinctionRate,
                double rho,
                unsigned numLociToSim,
                double geneBirthRate,
                double geneDeathRate,
                double transferRate,
                unsigned indPerPop,
                unsigned popSize,
                double genTime,
                int ng,
                double og,
                double ts);
        ~Simulator();

        void    setSpeciesTree(SpeciesTree *st) { spTree = st; }
        bool    gsaBDSim();
        bool    bdsaBDSim();
        bool    moranSpeciesSim();
        bool    coalescentSim();
        bool    simSpeciesTree();
        bool    simMoranSpeciesTree();
        bool    simSpeciesLociTrees();
        bool    simThreeTree();
        bool    simLocusGeneTrees();
        bool    gsaCheckStop();
        void    initializeSim();
        void    processGSASim();
        void    prepGSATreeForReconstruction();
        void    processSpTreeSim();
        void    graftOutgroup(Tree *tr, double trDepth);
        std::string    printSpeciesTreeNewick();
        std::string    printLocusTreeNewick(int i);
        std::string    printGeneTreeNewick(int i, int j);
        std::string    printExtantGeneTreeNewick(int i, int j);
        std::set<double, std::greater<double> > getEpochs();
};


#endif /* Simulator_h */
