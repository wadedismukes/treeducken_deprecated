//
//  GeneTree.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef GeneTree_h
#define GeneTree_h

#include "LocusTree.h"

class GeneTree : public Tree {
    private:
        unsigned numTips;
        unsigned individualsPerPop;
        unsigned popSize;
        double   generationTime;
    
    public:
                    GeneTree(MbRandom *rando, unsigned nt, unsigned ipp, unsigned ne, double genTime);
                    ~GeneTree();
        double      getTimeToNextEvent(int n); // what do you need to determine this?
        Node*       coalescentEvent(double t, Node *p, Node *q);
        bool        censorCoalescentProcess(double startTime, double stopTime, int contempSpIndx, int newSpIndx);
        void        initializeTree(std::unordered_set<int> extantLociIndx, std::multimap<int,double> locusDeathTimes);
        std::multimap<int,double> rescaleTimes(std::multimap<int, double> timeMap);
        void        rootCoalescentProcess(double startTime);
        void        recursiveRescaleTimes(Node *r, double add);
        void        setBranchLengths();
        void        setIndicesBySpecies(std::map<int,int> spToLocusMap);
        std::string printNewickTree();
        void        recGetNewickTree(Node *r, std::stringstream &ss);
        void        setTreeTipNames();
        void        addExtinctSpecies(double bt, int indx);
    
};

#endif /* GeneTree_hpp */
