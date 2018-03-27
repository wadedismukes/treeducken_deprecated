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
//TODO: fix something

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
        (*p)->setFlag(1);
    }
    setSampleFromFlags();
}


void Tree::setSampleFromFlags(){
    
    
    int flag;
    Node *q = NULL;
    for(std::vector<Node *>::iterator p=nodes.begin(); p!=nodes.end(); p++){
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

