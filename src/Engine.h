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
            double                      treeLength, treeNess, aveTipLen, treeDepth;
            
        public:
                                        TreeInfo(int idx) : treeIndx(idx), treeLength(0.0), treeNess(0.0), aveTipLen(0.0) {};
            std::string                 getWholeSpeciesTree();
            std::string                 getLocusTreeByIndx(int idx);
            std::string                 getGeneTreeByIndx(int idx);
            double                      getTreeLength() {return treeLength; }
            double                      getTreeNess() {return treeNess; }
            double                      getAveTipLen() {return aveTipLen; }
            double                      getTreeDepth() {return treeDepth; }
            
            
            void                        setWholeTreeStringInfo(std::string ts ) { speciesTree = ts; }
            void                        setLocusTreeByIndx(int indx, std::string ts) { locusTrees[indx] = ts; }
            void                        setGeneTreeByIndx(int indx, std::string ts) { geneTrees[indx] = ts; }
            void                        setTreeLength(double b) { treeLength = b; }
            void                        setTreeNess(double b) { treeNess = b; }
            void                        setAveTipLen(double b) {aveTipLen = b; }
            void                        setTreeDepth(double b) { treeDepth = b; }
            
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
        
                                Engine(std::string of, int mt, double sbr, double sdr, double gbr, double gdr, double lgtr, int ipp, int popsize, double genTime);
                                ~Engine();
        void                    doRunRun();
        void                    writeTreeFiles();
        TreeInfo                *findTreeByIndx(int i);
        
        
        
};

#endif /* Engine_h */
