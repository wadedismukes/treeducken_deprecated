//
//  Tree.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/7/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Tree_h
#define Tree_h

#include <string>
#include <vector>
#include "MbRandom.h"
#include <iostream>

class Node
{
    private:
        Node    *ldes;
        Node    *rdes;
        Node    *anc;
        Node    *sib;
        int     indx, Lindx;
        int     flag;
        std::string name;
        bool    isRoot;
        bool    isTip;
        bool    isExtant, isExtinct;
        bool    isDuplication;
        double  birthTime, deathTime;
        double  branchLength;
        
    public:
                Node();
                ~Node();
        void    setAsRoot(bool t) {isRoot = t; }
        void    setBirthTime(double bt) {birthTime = bt; }
        void    setIndx(unsigned i) {indx = i; }
        void    setIsTip(bool t) {isTip = t; }
        void    setDeathTime(double dt) {deathTime = dt; }
        void    setIsExtant(bool t) {isExtant = t; }
        void    setIsExtinct(bool t) {isExtinct = t; }
        void    setLdes(Node *l) {ldes = l; }
        void    setRdes(Node *r) {rdes = r; }
        void    setAnc(Node *a) {anc = a; }
        void    setSib(Node *s) {sib = s; }
        void    setName(std::string f) { name = f; }
        void    setBranchLength(double bl) {branchLength = bl; } 
        void    setFlag(int d) { flag = d; }
        void    setIndx(int i) {indx = i; }
        void    setLindx(int li ) {Lindx = li; }
        void    setIsDuplication(bool t) { isDuplication = t; }
    
        int     getFlag() {return flag; }
        Node*   getLdes() {return ldes; }
        Node*   getRdes() {return rdes; }
        Node*   getAnc() {return anc; }
        Node*   getSib() {return sib; }
        bool    getIsRoot() {return isRoot; }
        bool    getIsTip() {return isTip; }
        bool    getIsExtinct() {return isExtinct; }
        bool    getIsExtant() { return isExtant; }
        std::string getName() { return name; }
        double  getBranchLength() { return branchLength; }
        double  getDeathTime() {return deathTime; }
        double  getBirthTime() { return birthTime; }
        int     getIndex() {return indx; }
        int     getLindx() { return Lindx; }
        bool    getIsDuplication() { return isDuplication; }
};



class Tree
{
    protected:
        Node    *root;
        Node    *outgrp;
        std::vector<Node*> nodes;
        std::vector<Node*> extantNodes;
        unsigned numTaxa, numNodes, numTotalTips;
        unsigned numExtant, numExtinct;
        double  currentTime;
        MbRandom *rando;

    public:
                    Tree(MbRandom *p, unsigned numExtant, double cTime);
                    Tree(MbRandom *p, unsigned numTaxa);
                    ~Tree();
        void        setOutgroup(Node *og) { outgrp = og; }
        Node*       getOutgroup() { return outgrp; }
        Node*       getRoot() {return root; }
        void        setRoot(Node *r) { root = r; }
        double      getNumExtant() {return numExtant; }
        double      getNumExtinct() {return numExtinct; }
        int         getNodesSize() { return (int) nodes.size(); }
        double      getTotalTreeLength();
        double      getTreeDepth();
        double      getCurrentTime() {return currentTime; }
        double      getEndTime();
        void        rescaleTreeByOutgroupFrac(double outgroupFrac, double getTreeDepth);
        void        zeroAllFlags();
        void        setWholeTreeFlags();
        void        setExtantTreeFlags();
        void        setSampleFromFlags();
        void        getRootFromFlags();
        void        getExtantTree();
        void        setNewRootInfo(Node *newRoot, Node *outgroup, Node *oldRoot, double t);
        std::vector<Node*> getNodes() { return nodes; }
        std::vector<Node*> getExtantNodes() { return extantNodes; }
        
        void        reconstructTreeFromSim(Node *oRoot);
        void        reconstructLineageFromSim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter);

        virtual double  getTimeToNextEvent() { return 0.0; }
        virtual void    lineageBirthEvent() { return; }
        virtual void    lineageDeathEvent() { return; }
        virtual void    setTreeTipNames()  { return; }
        virtual void    ermEvent(double ct) { return; }
        virtual void    setBranchLengths() { return; }
        virtual std::string    printNewickTree() { return "t";}
        friend class Node;
        
};
#endif /* Tree_hpp */
