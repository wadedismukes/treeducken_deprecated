#include <utility>

#include <utility>

#ifndef Engine_h
#define Engine_h

#include "Simulator.h"
#include <iostream>
#include <fstream>
#include <regex>

/**
 * @brief Class for handling the trees and data about trees from the simulation
 *        run. Holds vectors of newick trees as well as various tree statistics.
 *
 */

class TreeInfo{
        private:
            std::string                 speciesTree;
            std::string                 extSpeciesTree;
            std::vector<std::string>    gsaTrees;
            std::vector<std::string>    locusTrees;
            std::vector<std::vector<std::string> >   geneTrees;
            std::vector<std::vector<std::string> >   extGeneTrees;
            double                      spTreeLength, spTreeNess, spAveTipLen, spTreeDepth;
            double                      extSpTreeLength, extSpTreeDepth;
            double                      loTreeLength, loTreeNess, loAveTipLen, loTreeDepth;
            double                      aveTMRCAGeneTree;
            int                         numTransfers, numDuplications{}, numLosses{};
            std::vector<double>         numGenerations;
    

    
        public:
                                        TreeInfo(int idx, int nl);
                                        ~TreeInfo(); 
            std::string                 getWholeSpeciesTree() { return speciesTree; }
            std::string                 getExtantSpeciesTree() { return extSpeciesTree; }
            std::string                 getLocusTreeByIndx(int idx) { return locusTrees[idx]; }
            std::string                 getGeneTreeByIndx(int Lidx, int idx) { return geneTrees[Lidx][idx]; }
            std::string                 getExtGeneTreeByIndx(int Lidx, int idx) { return extGeneTrees[Lidx][idx]; }
            double                      getSpeciesTreeLength() {return spTreeLength; }
            double                      getSpeciesTreeNess() {return spTreeNess; }
            double                      getSpeciesAveTipLen() {return spAveTipLen; }
            double                      getSpeciesTreeDepth() {return spTreeDepth; }
            double                      getLocusTreeLength() {return loTreeLength; }
            double                      getLocusTreeNess() {return loTreeNess; }
            double                      getLocusAveTipLen() {return loAveTipLen; }
            double                      getLocusTreeDepth() {return loTreeDepth; }
            double                      getAveTMRCAGeneTree() {return aveTMRCAGeneTree; }
            int                         getNumberTransfers() { return numTransfers; }
            int                         getNumberDuplications() { return numDuplications; }
            int                         getNumberLosses() { return numLosses; }
            std::vector<double>         getNumberGenerations() { return numGenerations; }
            double                      getNumberGenerationsByLindx(int lindx) { return numGenerations[lindx]; }
            double                      getExtSpeciesTreeDepth() { return extSpTreeDepth; }
        
            void                        setNumberTransfers(int d) { numTransfers = d; }
            void                        setNumberDuplications(int d) { numDuplications = d; }
            void                        setNumberLosses(int d) { numLosses = d; }
            void                        setNumberGenerations(std::vector<double> ng) { numGenerations = std::move(ng); }
            void                        setWholeTreeStringInfo(std::string ts ) { speciesTree = std::move(ts); }
            void                        setExtTreeStringInfo(std::string ts) { extSpeciesTree = std::move(ts); }
            void                        setLocusTreeByIndx(int indx, const std::string& ts) { locusTrees.push_back(ts); }
            void                        setGeneTreeByIndx(int Lindx, int indx, const std::string& ts) { geneTrees[Lindx].push_back(ts); }
            void                        setExtantGeneTreeByIndx(int Lindx, int indx, const std::string& ts) { extGeneTrees[Lindx].push_back(ts); }
            void                        setSpeciesTreeLength(double b) { spTreeLength = b; }
            void                        setSpeciesTreeNess(double b) { spTreeNess = b; }
            void                        setSpeciesAveTipLen(double b) {spAveTipLen = b; }
            void                        setSpeciesTreeDepth(double b) { spTreeDepth = b; }
            void                        setLocusTreeLength(double b) { loTreeLength = b; }
            void                        setLocusTreeNess(double b) { loTreeNess = b; }
            void                        setLocusAveTipLen(double b) {loAveTipLen = b; }
            void                        setLocusTreeDepth(double b) { loTreeDepth = b; }
            void                        setAveTMRCAGeneTree(double b) { aveTMRCAGeneTree = b; }
            void                        setExtSpeciesTreeDepth(double b) { extSpTreeDepth = b; }
            void                        writeTreeStatsFile(int spIndx, std::string ofp);      
            void                        writeExtantTreeFileInfo(int spIndx, std::string ofp);                
            void                        writeWholeTreeFileInfo(int spIndx, std::string ofp);
            void                        writeLocusTreeFileInfoByIndx(int spIndx, int indx, std::string ofp);
            void                        writeGeneTreeFileInfoByIndx(int spIndx, int Lindx, int indx, std::string ofp);
            void                        writeGeneTreeFileInfo(int spIndx, int Lindx, int numgene, std::string ofp);
            void                        writeExtGeneTreeFileInfo(int spIndx, int Lindx, int numgene, std::string ofp);
};      


/**
 * @brief The engine class captures the settings and uses those to run the
 *        simulation functions. The member functions of this class are wrappers
 *        for the simulation functions in the Simulator class.
 */

class Engine{
    
    private:
        
        std::string outfilename;
        std::string inputSpTree;
        std::vector<TreeInfo*> simSpeciesTrees;
        int                    simType;
        int                    numTaxa;
        int                    numSpeciesTrees;
        int                    numLoci, numGenes;
        double                 treescale;
        bool                   doScaleTree;
        double                 outgroupFrac;
        MbRandom               rando;
        // species tree paramters
        double                 spBirthRate, spDeathRate;
        double                 proportionToSample;
        // locus tree parameters
        double                 geneBirthRate, geneDeathRate, transferRate;
        // gene tree paramters
        int                    individidualsPerPop, populationSize;
        double                 generationTime;
        bool                   printOutputToScreen;
        
    public:
        
                                Engine(std::string of, int mt, double sbr, double sdr, double gbr, double gdr,
                                       double lgtr, int ipp,
                                       int popsize, double genTime, int sd1, int sd2, double ts, int reps, int ntax,
                                       int nloci, int ngen,
                                       double og, bool sout);
                                ~Engine();
        static unsigned int            countNewickLeaves(const std::string& stNewick);
        static std::string             stripCommentsFromNewickTree(const std::string& stNewick);
        static std::string             formatTipNamesFromNewickTree(const std::string& stNewick);
        void                    setInputSpeciesTree(const std::string& stNewick);
        std::string             getInputSpeciesTree() { return inputSpTree; }
        void                    doRunRun();
        void                    doRunSpTreeSet();
        void                    writeTreeFiles();
        TreeInfo                *findTreeByIndx(int i);
        void                    calcAverageRootAgeSpeciesTrees();
        SpeciesTree*            buildTreeFromNewick(const std::string& spTree);
        
};

#endif /* Engine_h */
