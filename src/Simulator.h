//
//  Simulator.hpp
//  treeducken
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
        bool        printSOUT;
        std::vector<SpeciesTree*>   gsaTrees;
        SpeciesTree*    spTree;
        LocusTree*      lociTree;
        std::vector<LocusTree*> locusTrees;
        GeneTree*       geneTree;
        std::vector<std::vector<GeneTree*> > geneTrees;

    public:
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
                double ts,
                bool sout);
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
        bool    moranCheckStop();

    void    processGSASim();
        void    prepGSATreeForReconstruction();
        void    processSpTreeSim();
        void    graftOutgroup(Tree *tr, double trDepth);
        double  calcSpeciesTreeDepth(); 
        double  calcExtantSpeciesTreeDepth();
        double  calcLocusTreeDepth(int i);
        int     findNumberTransfers();
        int     findNumberDuplications();
        int     findNumberLosses();
        std::vector<double>  findAveNumberGenerations();
        std::string    printSpeciesTreeNewick();
        std::string    printExtSpeciesTreeNewick();
        std::string    printLocusTreeNewick(int i);
        std::string    printGeneTreeNewick(int i, int j);
        std::string    printExtantGeneTreeNewick(int i, int j);
        std::set<double, std::greater<double> > getEpochs();
};


#endif /* Simulator_h */
