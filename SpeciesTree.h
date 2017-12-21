//
//  SpeciesTree.h
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/8/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef SpeciesTree_h
#define SpeciesTree_h
#include "MbRandom.h"
#include "Tree.h"
#include <sstream>
#include <map>

class SpeciesTree : public Tree
{
    private:
    
        double        speciationRate, extinctionRate;
        unsigned        extantStop;

    public:
                      SpeciesTree(MbRandom *p, unsigned numTaxa, double curTime, double specRate, double extRate);
                      SpeciesTree(MbRandom *p, unsigned numTaxa);
                      ~SpeciesTree();
        void          setSpeciationRate(double sr) {speciationRate = sr; }
        void          setExtinctionRate(double er) {extinctionRate = er; }

        // tree-building functions
        double        getTimeToNextEvent();
        void          lineageBirthEvent(unsigned indx);
        void          lineageDeathEvent(unsigned indx);
        void          ermEvent(double curTime);
        void          setNewLineageInfo(unsigned indx, Node *r, Node *l);
    
        // set node parameters across tree
        void          setBranchLengths();
        void          setPresentTime(double currentT);
        void          setTreeTipNames();
        void          recTipNamer(Node *p, unsigned &extinctCount, unsigned &tipCount);
    
        std::string   printNewickTree();
        void          recGetNewickTree(Node *p, std::stringstream &ss);
    
        // simulation functions
        void          setGSATipTreeFlags();
        void          reconstructTreeFromGSASim(Node *oRoot);
        void          popNodes();
        void          recPopNodes(Node *p);
        void          reconstructLineageFromGSASim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter);
        void          setSampleFromFlags();
        std::map<int,double>        getBirthTimesFromNodes();
        std::map<int,double>        getDeathTimesFromNodes();
        double                      getCurrentTimeFromExtant() {return extantNodes[0]->getDeathTime();}
        bool                        getIsExtantFromIndx(int indx) { return nodes[indx]->getIsExtant(); }
        bool                        macroEvent(int indx);
        std::pair<int, int>         preorderTraversalStep(int index);
};


#endif /* SpeciesTree_h */
