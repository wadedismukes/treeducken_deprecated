//
//  Engine.h
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 1/28/18.
//  Copyright © 2018 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Engine_hpp
#define Engine_hpp
#include "Simulator.h"
#include <iostream>
#include <fstream>
/*
 tree info class
                */

class TreeInfo{
        private:
            int                         treeIndx;
            std::string                 speciesTree;
            std::vector<std::string>    gsaTrees;
            std::vector<std::string>    locusTrees;
            std::vector<std::string>    geneTrees;
            double                      spTreeLength, spTreeNess, spAveTipLen, spTreeDepth;
            double                      loTreeLength, loTreeNess, loAveTipLen, loTreeDepth;
            double                      aveTMRCAGeneTree;

    
        public:
                                        TreeInfo(int idx) : treeIndx(idx), spTreeLength(0.0), spTreeNess(0.0), spAveTipLen(0.0), loTreeLength(0.0), loTreeNess(0.0), loAveTipLen(0.0) {};
            std::string                 getWholeSpeciesTree() {return speciesTree; }
            std::string                 getLocusTreeByIndx(int idx) { return locusTrees[idx]; }
            std::string                 getGeneTreeByIndx(int idx) { return geneTrees[idx]; }
            double                      getSpeciesTreeLength() {return spTreeLength; }
            double                      getSpeciesTreeNess() {return spTreeNess; }
            double                      getSpeciesAveTipLen() {return spAveTipLen; }
            double                      getSpeciesTreeDepth() {return spTreeDepth; }
            double                      getLocusTreeLength() {return loTreeLength; }
            double                      getLocusTreeNess() {return loTreeNess; }
            double                      getLocusAveTipLen() {return loAveTipLen; }
            double                      getLocusTreeDepth() {return loTreeDepth; }
            double                      getAveTMRCAGeneTree() {return aveTMRCAGeneTree; }
    
        
    
            void                        setWholeTreeStringInfo(std::string ts ) { speciesTree = ts; }
            void                        setLocusTreeByIndx(int indx, std::string ts) { locusTrees[indx] = ts; }
            void                        setGeneTreeByIndx(int indx, std::string ts) { geneTrees[indx] = ts; }
            void                        setSpeciesTreeLength(double b) { spTreeLength = b; }
            void                        setSpeciesTreeNess(double b) { spTreeNess = b; }
            void                        setSpeciesAveTipLen(double b) {spAveTipLen = b; }
            void                        setSpeciesTreeDepth(double b) { spTreeDepth = b; }
            void                        setLocusTreeLength(double b) { loTreeLength = b; }
            void                        setLocusTreeNess(double b) { loTreeNess = b; }
            void                        setLocusAveTipLen(double b) {loAveTipLen = b; }
            void                        setLocusTreeDepth(double b) { loTreeDepth = b; }
            void                        setAveTMRCAGeneTree(double b) {aveTMRCAGeneTree = b; }
    
    
            void                        writeWholeTreeFileInfo(std::string ofp);
            void                        writeLocusTreeFileInfoByIndx(int indx, std::string ofp);
            void                        writeGeneTreeFileInfoByIndx(int indx, std::string ofp);

};


/*
  engine class
                */

class Engine{
    
    private:
        
        std::string outfilename;
        std::vector<TreeInfo*> simSpeciesTrees;
        int                    simType;
        int                    seedset;
        int                    numTaxa;
        int                    numSpeciesTrees;
        int                    numLoci, numGenes;
        double                 treescale;
        bool                   doScaleTree;
    
        MbRandom               rando;
        // species tree paramters
        double                 spBirthRate, spDeathRate;
        double                 proportionToSample;
        // locus tree parameters
        double                 geneBirthRate, geneDeathRate, transferRate;
        // gene tree paramters
        int                    individidualsPerPop, populationSize;
        double                 generationTime;
        
    public:
        
                                Engine(std::string of,
                                       int mt,
                                       double sbr,
                                       double sdr,
                                       double gbr,
                                       double gdr,
                                       double lgtr,
                                       int ipp,
                                       int popsize,
                                       double genTime,
                                       int sd1,
                                       int sd2,
                                       double treescale,
                                       int reps,
                                       int numTaxa,
                                       int nloci);
                                ~Engine();
        void                    doRunRun();
        void                    writeTreeFiles();
        TreeInfo                *findTreeByIndx(int i);
        void                    calcAverageRootAgeSpeciesTrees();
    
        
        
};

#endif /* Engine_h */
