//
//  SpeciesTree.h
//  treeducken
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
#include <set>

class SpeciesTree : public Tree
{
    private:
    
        double        speciationRate{}, extinctionRate{};
        unsigned      extantStop;

    public:
                      SpeciesTree(MbRandom *p, unsigned numTaxa, double br, double dr);
                      SpeciesTree(MbRandom *p, unsigned numTaxa);
        virtual       ~SpeciesTree();
        void          setSpeciationRate(double sr) {speciationRate = sr; }
        void          setExtinctionRate(double er) {extinctionRate = er; }

        // tree-building functions
        virtual double        getTimeToNextEvent(); 
        double        getTimeToNextEventMoran();
        virtual void          lineageBirthEvent(unsigned indx);
        virtual void          lineageDeathEvent(unsigned indx);
        void          ermEvent(double curTime);
        void          moranEvent(double curTime); 
        void          setNewLineageInfo(unsigned indx, Node *r, Node *l);
        void          initializeMoranProcess(unsigned numTaxa);
    
        // set node parameters across tree
        void          setBranchLengths();
        void          setPresentTime(double currentT);
        void          setTreeTipNames();
        static void          recTipNamer(Node *p, unsigned &extinctCount, unsigned &tipCount);
        

        std::string   printNewickTree();
        std::string   printExtNewickTree();
        static void          recGetNewickTree(Node *p, std::stringstream &ss);
    
        // simulation functions
        void          setGSATipTreeFlags();
        void          reconstructTreeFromGSASim(Node *oRoot);
        void          setTreeInfo();
        void          popNodes();
        void          recPopNodes(Node *p);
        void          reconstructLineageFromGSASim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter);
      //  void          setSampleFromFlags();
        void          initializeMoranProcess(); // TODO: Write this function
        std::map<int,double>        getBirthTimesFromNodes();
        std::map<int,double>        getDeathTimesFromNodes();
        double                      getCurrentTimeFromExtant() {return extantNodes[0]->getDeathTime();}
        bool                        getIsExtantFromIndx(int indx) { return nodes[indx]->getIsExtant(); }
        bool                        macroEvent(int indx);
        std::pair<int, int>         tipwiseStep(int index);
        int                         rootwiseStep(int index);
};


#endif /* SpeciesTree_h */
