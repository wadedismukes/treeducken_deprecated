//
//  LocusTree.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/13/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef LocusTree_h
#define LocusTree_h
#include "SpeciesTree.h"
#include <algorithm>

class LocusTree : public Tree 
{
    private:
        double geneBirthRate, geneDeathRate, transferRate;
        double currentTime;
        double stopTime;
        unsigned numTaxa;

    public:
        LocusTree(MbRandom *rando, unsigned nt, double stop, double gbr, double gdr, double lgtr);
        ~LocusTree();
        double  getTimeToNextEvent();
        void    lineageBirthEvent(int indx);
        void    lineageDeathEvent(int indx);
        void    setNewLineageInfo(int indx, Node *r, Node *s);
        void    lineageTransferEvent(int indx);
        void    ermEvent(double ct);
    
    
        int     speciationEvent(int indx, double time);
        void    extinctionEvent(int indx, double time);
        void    setNewIndices(int indx, std::pair<int,int> sibs, int count);
        std::string   printNewickTree();
        void    setTreeTipNames();
        void    recTipNamer(Node *p, unsigned &copyNumber);
        void    recGetNewickTree(Node *r, std::stringstream &ss);
        void    setBranchLengths();
        void    setPresentTime(double currentT);
        void    setStopTime(double st) {stopTime = st; }
    

};
#endif /* LocusTree_hpp */
