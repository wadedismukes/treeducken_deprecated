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

void GeneTree::initializeTree(std::vector< std::vector<int> > extantLociInd, double presentTime){
    if(!(nodes.empty())){
        nodes.clear();
    }
    if(!(extantNodes.empty())){
        extantNodes.clear();
    }
    Node *p;
    int numLociInPresnt = extantLociInd[0].size();
    for(int i = 0; i < numLociInPresnt; i++){
        for(int j = 0; j < individualsPerPop; j++){
            p = new Node();
            p->setDeathTime(presentTime);
            p->setLindx(extantLociInd[0][i]);
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
    ct = -log(rando->uniformRv()) / (lambda);
 //   std::cout << "coal time " << ct << std::endl;
    return ct;
}

bool GeneTree::censorCoalescentProcess(double startTime, double stopTime, int contempSpeciesIndx, int ancSpIndx, bool chck){
    int leftInd, rightInd;
    int leftIndExtN, rightIndExtN;
    int extIndx;
    Node *l, *r;
    Node *n;
    double t = startTime;
    bool allCoalesced = false;
    // search extantNodes for members with Lindx = contempSpecisIndx 
    std::vector<int> indInExtNodes;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        if((*it)->getLindx() == contempSpeciesIndx){
            extIndx = std::distance(extantNodes.begin(), it);
            indInExtNodes.push_back(extIndx);
        }
    }
    std::cout << contempSpeciesIndx << "   ($)$)%    " << indInExtNodes.size() << std::endl;
    if(indInExtNodes.size() > 1){
        while(t > stopTime){
            t -= getTimeToNextEvent(indInExtNodes.size());
            std::cout << t << std::endl;
            if(t < stopTime){
                if(chck){
                    t = stopTime;
                    allCoalesced = true;
                    break;
                }
                else{
                    t = stopTime;
                    allCoalesced = false;
                    break;
                }
            }

            rightInd = rando->uniformRv(0, indInExtNodes.size() - 1);
            rightIndExtN = indInExtNodes[rightInd];
            r = extantNodes[rightIndExtN];
            indInExtNodes.erase(indInExtNodes.begin() + rightInd);

            leftInd = rando->uniformRv(0, indInExtNodes.size() - 1);
            leftIndExtN = indInExtNodes[leftInd];
            l = extantNodes[leftIndExtN];
            indInExtNodes.erase(indInExtNodes.begin() + leftInd);

            n = coalescentEvent(t, l, r);

            iter_swap(extantNodes.begin() + rightIndExtN, extantNodes.end() - 1);
            iter_swap(extantNodes.begin() + leftIndExtN, extantNodes.end() - 2);
            extantNodes.insert(extantNodes.begin(), n);
            if(!(indInExtNodes.empty())){
                for(int p = 0; p < indInExtNodes.size(); p++)
                    indInExtNodes[p] += 1;
            }
            indInExtNodes.push_back(0);
            
            extantNodes.erase(extantNodes.end() - 2, extantNodes.end());
            std::cout << "size of extantNodes " << extantNodes.size() << std::endl;
            if(indInExtNodes.size() == 1){
                allCoalesced = true;
                break;
            }
        }
    }
    else if (indInExtNodes.size() == 1){
        t = stopTime;
        allCoalesced = true;
        extantNodes[indInExtNodes[0]]->setLindx(ancSpIndx);
    }
    else{
        // std::cout << "this shouldn't even be happnening" << std::endl;
        allCoalesced = true;
    }
    
    if(allCoalesced == true){
        for(int i = 0; i < indInExtNodes.size(); ++i){
            std::cout << "**************" << std::endl;
            extantNodes[indInExtNodes[i]]->setLindx(ancSpIndx);
        }
    }


    return allCoalesced;
}

Node* GeneTree::coalescentEvent(double t, Node *p, Node *q){
    Node *n = new Node();
    n->setDeathTime(t);
    n->setLdes(p);
    n->setRdes(q);
    n->setIsExtant(false);
    n->setIsTip(false);
    n->setIsExtinct(false);
    n->setLindx(p->getLindx());
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
        p.second = ((*it).second);
        rescaledTimeMap.insert(p);
    }
    
    return rescaledTimeMap;
    
}

void GeneTree::rootCoalescentProcess(double startTime){
    int leftInd, rightInd;
    int leftIndExtN, rightIndExtN;
    int extIndx;
    Node *l, *r;
    Node *n;
    double t = startTime;
    // search extantNodes for members with Lindx = contempSpecisIndx 
    std::vector<int> indInExtNodes;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setLindx(0);
    }
    while(extantNodes.size() > 1){
        t -= getTimeToNextEvent(extantNodes.size());

        rightInd = rando->uniformRv(0, extantNodes.size() - 1);
        r = extantNodes[rightInd];
        extantNodes.erase(extantNodes.begin() + rightInd);

        leftInd = rando->uniformRv(0, extantNodes.size() - 1);
        l = extantNodes[leftInd];
        extantNodes.erase(extantNodes.begin() + leftInd);

        n = coalescentEvent(t, l, r);
        extantNodes.push_back(n);
    }


   

    this->setRoot(extantNodes[0]);
//    if(extantNodes[0]->getBirthTime() < 0){
//        time = std::abs(this->getRoot()->getDeathTime() * 2);
//        this->getRoot()->setBirthTime(0.0);
//        this->getRoot()->setDeathTime(this->getRoot()->getDeathTime() + time);
//        this->recursiveRescaleTimes(this->getRoot(), time);
//    }
    this->setBranchLengths();
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
        p->setLindx(indx);
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
        indx = (*it)->getLindx();
        spIndx = spToLocusMap.find(indx)->second;
        (*it)->setIndx(spIndx);
    }
    this->setTreeTipNames();
}

std::string GeneTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string geneTreeString = ss.str();
    return geneTreeString;
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
    int indNumber = 0;
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
            indNumber++;
            tn << indNumber;
            name += "_" + tn.str();
            tn.clear();
            tn.str(std::string());
            (*it)->setName(name);
            if(indNumber == individualsPerPop)
                indNumber = 0;
        }

    }

}
