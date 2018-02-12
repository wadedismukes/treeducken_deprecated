//
//  GeneTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "GeneTree.h"
#include <iostream>
GeneTree::GeneTree(MbRandom *p, unsigned nt, unsigned ipp, unsigned ne, double genTime) : Tree(p, nt){
    rando = p;
    numTaxa = nt;
    individualsPerPop = ipp;
    popSize = ne;
    generationTime = genTime;
    
}

GeneTree::~GeneTree(){
    
}

void GeneTree::initializeTree(std::unordered_set<int> extantLociInd, std::multimap<int,double> locusDeathMap){
    if(!(nodes.empty())){
        nodes.clear();
    }
    if(!(extantNodes.empty())){
        extantNodes.clear();
    }
    Node *p;
    for(std::unordered_set<int>::iterator it = extantLociInd.begin(); it != extantLociInd.end(); ++it){
        for(int j = 0; j < individualsPerPop; j++){
            p = new Node();
            p->setDeathTime(locusDeathMap.find(*it)->second);
            p->setIndx(*it);
            p->setLdes(NULL);
            p->setRdes(NULL);
            p->setAnc(NULL);
            p->setIsExtant(true);
            p->setIsTip(true);
            p->setIsExtinct(false);
            extantNodes.push_back(p);
            nodes.push_back(p);
        }
    }
}

double GeneTree::getTimeToNextEvent(int n){
    double ct;
    double lambda = (double)(n * (n - 1)) / (2 * popSize) ;
    ct = -log(rando->uniformRv()) * (lambda);
    
    return ct;
}


bool GeneTree::censorCoalescentProcess(double startTime, double stopTime, int contempSpeciesIndx, int ancSpIndx){
    int indx;
    Node *p;
    std::vector<Node*> coalescingNodes;
    double time = startTime;
    int leftInd, rightInd;
    bool reachedTime = false;
    Node *r, *l;
    std::vector<Node*>::iterator it = extantNodes.begin();
    for(; it != extantNodes.end(); ){
        indx = (*it)->getIndex();
        if(indx == contempSpeciesIndx){
            p = *it;
            coalescingNodes.push_back(p);
            it = extantNodes.erase(it);

        }
        else{
            ++it;
        }
    
    }
    if(coalescingNodes.size() < 1)
        return reachedTime = true;
    else{
        it = coalescingNodes.begin();
    }
    
    while(coalescingNodes.size() > 1){
        time -= getTimeToNextEvent((int)coalescingNodes.size());

        if(time < stopTime){
            reachedTime = true;
            break;
        }
        rightInd = rando->uniformRv(0, coalescingNodes.size() - 1);
        std::advance(it, rightInd);
        r = *it;
        coalescingNodes.erase(it);
        it = coalescingNodes.begin();

        leftInd = rando->uniformRv(0, coalescingNodes.size() - 1);
        std::advance(it, leftInd);
        l = *it;
        coalescingNodes.erase(it);
        it = coalescingNodes.begin();
        
        p = coalescentEvent(time, l, r);
        coalescingNodes.push_back(p);


    }
    
    it = coalescingNodes.begin();
    if(!(reachedTime)){
        for(; it != coalescingNodes.end(); ++it){
            time = stopTime;
            (*it)->setIndx(ancSpIndx);
            extantNodes.push_back((*it));
        }
    }
    else{
        for(; it != coalescingNodes.end(); ++it){
            (*it)->setIndx(contempSpeciesIndx);
            extantNodes.push_back((*it));
        }

    }
    coalescingNodes.clear();
    return reachedTime;
}

Node* GeneTree::coalescentEvent(double t, Node *p, Node *q){
    Node *n = new Node();
    n->setDeathTime(t);
    n->setLdes(p);
    n->setRdes(q);
    n->setIsExtant(false);
    n->setIsTip(false);
    n->setIsExtinct(false);
    n->setIndx(p->getIndex());
    nodes.push_back(n);
    
    p->setBirthTime(t);
    p->setAnc(n);
    p->setSib(q);
    
    q->setBirthTime(t);
    q->setAnc(n);
    q->setSib(p);
    
    
    return n;
}

std::multimap<int, double> GeneTree::rescaleTimes(std::multimap<int, double> timeMap){
    std::multimap<int, double> rescaledTimeMap;
    std::pair<int, double> p;
    for(std::multimap<int, double>::iterator it = timeMap.begin(); it != timeMap.end(); ++it){
        p.first = (*it).first;
        p.second = ((*it).second * generationTime) ;
        rescaledTimeMap.insert(p);
    }
    
    return rescaledTimeMap;
    
}

void GeneTree::rootCoalescentProcess(double startTime){
    Node *r, *l, *p;
    double time = startTime;
    int rightInd, leftInd;
    std::vector<Node*>::iterator it = extantNodes.begin();
    while(extantNodes.size() > 1){
        time -= getTimeToNextEvent((int)extantNodes.size());

        rightInd = rando->uniformRv(0, extantNodes.size() - 1);
        std::advance(it, rightInd);
        r = *it;
        extantNodes.erase(it);
        it = extantNodes.begin();
        
        leftInd = rando->uniformRv(0, extantNodes.size() - 1);
        std::advance(it, leftInd);
        l = *it;
        extantNodes.erase(it);
        it = extantNodes.begin();
        
        p = coalescentEvent(time, l, r);
        extantNodes.push_back(p);
    }
    
    
    this->setRoot(extantNodes[0]);
    extantNodes.clear();
//    if(extantNodes[0]->getBirthTime() < 0){
//        time = std::abs(this->getRoot()->getDeathTime() * 2);
//        this->getRoot()->setBirthTime(0.0);
//        this->getRoot()->setDeathTime(this->getRoot()->getDeathTime() + time);
//        this->recursiveRescaleTimes(this->getRoot(), time);
//    }
    this->setBranchLengths();
    this->setTreeTipNames();
}


void GeneTree::recursiveRescaleTimes(Node* r, double add){
    if(r != NULL){
        if( r->getRdes() == NULL){
            r->setBirthTime(r->getBirthTime() + add);
            r->setDeathTime(r->getDeathTime() + add);
        }
        else{

            r->getLdes()->setBirthTime(r->getLdes()->getBirthTime() + add);
            r->getLdes()->setDeathTime(r->getLdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getLdes(), add);

            r->getRdes()->setBirthTime(r->getRdes()->getBirthTime() + add);
            r->getRdes()->setDeathTime(r->getRdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getRdes(), add);

        }
    }
}

void GeneTree::setBranchLengths(){
    double brlen;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        brlen = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(brlen);
    }
}

void GeneTree::addExtinctSpecies(double bt, int indx){
    Node *p;
    for(int i = 0; i < individualsPerPop; i++){
        p = new Node();
        p->setDeathTime(bt);
        p->setIndx(indx);
        p->setLdes(NULL);
        p->setRdes(NULL);
        p->setAnc(NULL);
        p->setIsExtant(false);
        p->setIsTip(true);
        p->setIsExtinct(true);
        extantNodes.push_back(p);
        nodes.push_back(p);
    }
}


void GeneTree::setIndicesBySpecies(std::map<int, int> spToLocusMap){
    int indx;
    int spIndx;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        indx = (*it)->getIndex();
        spIndx = spToLocusMap.find(indx)->second;
        (*it)->setIndx(spIndx);
    }
}

std::string GeneTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string spTree = ss.str();
    return spTree;
}

void GeneTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL )
            ss << p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ")";
        }
    }
}

void GeneTree::setTreeTipNames(){
    std::stringstream tn;
    std::string name;
    int indx;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if((*it)->getIsTip()){
            indx  = (*it)->getIndex();
            tn << indx;
            name = tn.str();
            tn.clear();
            tn.str(std::string());

        }

    }

}