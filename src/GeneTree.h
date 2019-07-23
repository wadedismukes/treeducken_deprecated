#ifndef GeneTree_h
#define GeneTree_h

#include "LocusTree.h"
#include <algorithm>
/**
 * @brief GeneTree class which is a child of the Tree class. GeneTree is produced via functions within the
 *
 */
class GeneTree : public Tree {
    private:
        unsigned individualsPerPop; //! individuals per population to be sampled in coalescent functions
        unsigned popSize; //! population size used in the getCoalTime function
        double   generationTime; //! specified in generations per year

    public:
                    GeneTree(MbRandom *rando, unsigned nt, unsigned ipp, unsigned ne, double genTime);
        virtual     ~GeneTree();
        double      getCoalTime(int n);
        Node*       coalescentEvent(double t, Node *p, Node *q);
        bool        censorCoalescentProcess(double startTime, double stopTime, int contempSpIndx, int newSpIndx, bool chck);
        void        initializeTree(std::vector< std::vector<int> > extantLociIndx, double presentTime);
        std::multimap<int,double> rescaleTimes(std::multimap<int, double> timeMap);
        void        rootCoalescentProcess(double startTime, double ogf);
        void        recursiveRescaleTimes(Node *r, double add);
        void        setBranchLengths();
        void        setIndicesBySpecies(std::map<int,int> spToLocusMap);
        std::string printNewickTree();
        std::string printExtantNewickTree();
        void        recGetNewickTree(Node *r, std::stringstream &ss);
        void        recGetExtNewickTree(Node *r, std::stringstream &ss);
        void        setTreeTipNames();
        void        addExtinctSpecies(double bt, int indx);

};

#endif /* GeneTree_hpp */
