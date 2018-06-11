//
//  Tree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/7/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Tree.h"
#include <vector>
#include <string>

Node::Node()
{
    ldes = NULL;
    rdes = NULL;
    anc = NULL;
    sib = NULL;
    indx = -1;
    Lindx = -1;
    flag = -1;
    isRoot = false;
    isTip = false;
    isExtant = false;
    isDuplication = false;
    isExtinct = false;
    branchLength = 0.0;
    birthTime = 0.0;
    deathTime = 0.0;

    
}

Node::~Node(){
    
}




Tree::Tree(MbRandom *p, unsigned numExta, double curTime){
    rando = p;
    numNodes = 0;
    outgrp = nullptr;
    // intialize tree with root
    Node *r = new Node();
    r->setAsRoot(true);
    r->setBirthTime(0.0);
    r->setIndx(0);
    root = r;
    nodes.push_back(r);
    extantNodes.push_back(r);
    numExtant = 1;
    numTaxa = numExta;
    numExtinct = 0;
    currentTime = curTime;
    
}

Tree::Tree(MbRandom *p, unsigned numTax){
    numTaxa = numTax;
    rando = p;
    numNodes = 2 * numTax - 1;
    // intialize tree with root
    Node *r = new Node();
    r->setAsRoot(true);
    r->setBirthTime(0.0);
    r->setIndx(0);
    r->setIsExtant(true);
    root = r;
    nodes.push_back(r);
    extantNodes.push_back(r);
    numExtant = 1;
    currentTime = 0.0;
}

Tree::~Tree()
{
    if(!(nodes.empty() == false))
        nodes.clear();
    if(!(extantNodes.empty()))
        extantNodes.clear();
}


void Tree::zeroAllFlags(){
    for(std::vector<Node*>::iterator it=nodes.begin(); it!=nodes.end(); it++){
        (*it)->setFlag(0);
    }
}

void Tree::setWholeTreeFlags(){
    this->zeroAllFlags();
    numTotalTips = 0;
    for(std::vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); ++p){
        if((*p)->getIsTip()){
            (*p)->setFlag(1);
            numTotalTips++;
        }
    }
    setSampleFromFlags();
}


void Tree::setExtantTreeFlags(){
    this->zeroAllFlags();
    numTotalTips = 0;
    for(std::vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); ++p){
        if((*p)->getIsExtant())
            (*p)->setFlag(1);
    }
    this->setSampleFromFlags();
}


void Tree::setSampleFromFlags(){
    int flag;
    Node *q = NULL;
    for(std::vector<Node *>::iterator p=nodes.begin(); p!=nodes.end(); ++p){
        if((*p)->getIsTip()){
            flag = (*p)->getFlag();
            q = (*p);
            if(flag == 1){
                do{
                    q = q->getAnc();
                    flag = q->getFlag();
                    flag++;
                    q->setFlag(flag);
                }while (q->getIsRoot() == false && flag < 2);
            }
        }
    }
}


double Tree::getTotalTreeLength(){
    double sum = 0.0;
    for(std::vector<Node*>::iterator p = nodes.begin(); p != nodes.end(); ++p){
        Node *n = (*p);
        sum += n->getBranchLength();
    }
    return sum;
}

double Tree::getTreeDepth(){
    double td = 0.0;
    Node *r = this->getRoot();
    while(r->getIsTip() == false){
        td += r->getBranchLength();
        if(!(r->getLdes()->getIsExtinct()))
            r = r->getLdes();
        else
            r = r->getRdes();
    }
    return td;
}

void Tree::reconstructTreeFromSim(Node *oRoot){
    Node *n = new Node();
    unsigned tipCounter = numExtant;
    unsigned intNodeCounter = 0;
    reconstructLineageFromSim(n, oRoot, tipCounter, intNodeCounter);
    delete n;
}

void Tree::reconstructLineageFromSim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter){
    Node *p;
    bool rootN = prevN->getIsRoot();
    double brlen = prevN->getBranchLength();
    int oFlag = prevN->getFlag();
    if(prevN->getIsTip() && oFlag == 1){
        // need to recalculate branchlength
        Node *prevAnc = prevN->getAnc();
        int ancFlag = prevAnc->getFlag();
        if(ancFlag == 1){
            brlen += prevAnc->getBranchLength();
            while(!prevAnc->getIsRoot() && ancFlag < 2){
                prevAnc = prevAnc->getAnc();
                ancFlag = prevAnc->getFlag();
                if(ancFlag == 1)
                    brlen += prevAnc->getBranchLength();
            }
        }
        
        p = new Node();
        tipCounter++;
        p->setBranchLength(brlen);
        p->setIsTip(true);
        p->setName(prevN->getName());
        p->setBirthTime(prevN->getBirthTime());
        p->setDeathTime(prevN->getDeathTime());
        p->setIsExtant(prevN->getIsExtant());
        p->setIsExtinct(prevN->getIsExtinct());
        p->setAnc(currN);
        if(currN->getLdes() == NULL)
            currN->setLdes(p);
        else if(currN->getRdes() == NULL)
            currN->setRdes(p);
        else{
            std::cerr << "ERROR: Problem adding a tip to the tree!" << std::endl;
            exit(1);
        }
        
    }
    else{
        if(oFlag > 1){
            Node *s1 = new Node();
            intNodeCounter++;
            if(prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromSim(s1, prevN->getLdes(), tipCounter, intNodeCounter);
            if(prevN->getRdes()->getFlag() > 0)
                reconstructLineageFromSim(s1, prevN->getRdes(), tipCounter, intNodeCounter);
            
            
            if(rootN == false){
                Node *prevAnc = prevN->getAnc();
                int ancFlag = prevAnc->getFlag();
                if(ancFlag == 1){
                    brlen += prevAnc->getBranchLength();
                    while(!(prevAnc)->getIsRoot() && ancFlag < 2){
                        prevAnc = prevAnc->getAnc();
                        ancFlag = prevAnc->getFlag();
                        if(ancFlag == 1)
                            brlen += prevAnc->getBranchLength();
                    }
                }
                
                if(currN != NULL){
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                    s1->setAnc(currN);
                    if(currN->getLdes() == NULL)
                        currN->setLdes(s1);
                    else if(currN->getRdes() == NULL)
                        currN->setRdes(s1);
                    else{
                        std::cerr << "ERROR: Probem adding an internal node to the tree" << std::endl;
                        exit(1);
                    }
                }
                else{
                    s1->setAsRoot(true);
                    setRoot(s1);
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                }
                
            }
            else{
                s1->setAsRoot(true);
                setRoot(s1);
                s1->setBranchLength(0.0);
                s1->setBirthTime(prevN->getBirthTime());
                s1->setDeathTime(prevN->getDeathTime());
            }

        }
        else if(oFlag == 1){
            if(prevN->getRdes()->getFlag() == 0 && prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromSim(currN, prevN->getLdes(), tipCounter, intNodeCounter);
            else
                reconstructLineageFromSim(currN, prevN->getRdes(), tipCounter, intNodeCounter);
        }
    }
}
// Gene tree version only 
void Tree::getRootFromFlags(){
    Node *p;

	setExtantTreeFlags();
    int numNodes = nodes.size() - 1;
	for(int i=numNodes; i > 0; i--){
		p = nodes[i];
		if(p->getFlag() >= 2){
			root = p;
			root->setAsRoot(true);
			break;
		}
		
    }
}

void Tree::rescaleTreeByOutgroupFrac(double outgroupFrac, double treeDepth){
    double birthTime, deathTime;
    double rescaleFactor = outgroupFrac * treeDepth;
    for(std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); ++it){
        birthTime = (*it)->getBirthTime();
        deathTime = (*it)->getDeathTime();
        
        (*it)->setBirthTime(birthTime + rescaleFactor);
        (*it)->setDeathTime(deathTime + rescaleFactor);

        (*it)->setBranchLength((*it)->getDeathTime() - (*it)->getBirthTime());
    }
}

void Tree::setNewRootInfo(Node *rootN, Node *outgroupN, Node *currRoot, double t){
   // rootN->setBirthTime(0.0);
    rootN->setDeathTime(currRoot->getBirthTime());
    rootN->setBranchLength(rootN->getDeathTime() - rootN->getBirthTime());
    rootN->setAsRoot(true);
    rootN->setLdes(currRoot);
    rootN->setRdes(outgroupN);
    // nodes.push_back(rootN);
    this->setRoot(rootN);

    currRoot->setAsRoot(false);
    currRoot->setAnc(rootN);
    
    outgroupN->setName("OUT");
    outgroupN->setBirthTime(currRoot->getBirthTime());
    outgroupN->setDeathTime(t);
    outgroupN->setIsTip(true);
    outgroupN->setBranchLength(outgroupN->getDeathTime() - outgroupN->getBirthTime());
    outgroupN->setIsExtant(true);
    outgroupN->setAnc(rootN);
    outgroupN->setLdes(NULL);
    outgroupN->setRdes(NULL);
    this->setOutgroup(outgroupN);
    // nodes.push_back(outgroupN);
}

double Tree::getEndTime(){
    double tipDtime = 0.0;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsTip() && (*it)->getIsExtant()){
            tipDtime = (*it)->getDeathTime();
            break;
        }
    }

    return tipDtime;
}