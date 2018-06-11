//
//  SpeciesTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/8/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "SpeciesTree.h"
#include <iostream>

SpeciesTree::SpeciesTree(MbRandom *p, unsigned numTaxa, double ct, double br, double dr) : Tree(p, numTaxa, 0.0){
    currentTime = 0.0;
    rando = p;
    extantStop = numTaxa;
    speciationRate = br;
    extinctionRate = dr;
    
}

SpeciesTree::SpeciesTree(MbRandom *p, unsigned numTaxa) : Tree(p, numTaxa){
    rando = p;
    extantStop = numTaxa;
}

SpeciesTree::~SpeciesTree(){
    
}

double SpeciesTree::getTimeToNextEvent(){
    double sumrt = speciationRate + extinctionRate;
    double returnTime = 0.0;
    
    returnTime = -log(rando->uniformRv()) / (double(numExtant) * sumrt);
    return returnTime;
}

void SpeciesTree::lineageBirthEvent(unsigned indx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    setNewLineageInfo(indx, right, sis);
}

void SpeciesTree::lineageDeathEvent(unsigned int indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void SpeciesTree::ermEvent(double cTime){
    currentTime = cTime;
    int nodeInd = rando->discreteUniformRv(0, numExtant - 1);
    double relBr = speciationRate / (speciationRate + extinctionRate);
    bool isBirth = (rando->uniformRv() < relBr ? true : false);
    if(isBirth)
        lineageBirthEvent(nodeInd);
    else
        lineageDeathEvent(nodeInd);
}

void SpeciesTree::setNewLineageInfo(unsigned int indx, Node *r, Node *l){
    extantNodes[indx]->setLdes(l);
    extantNodes[indx]->setRdes(r);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsExtant(false);
    
    r->setLdes(NULL);
    r->setRdes(NULL);
    r->setSib(l);
    r->setAnc(extantNodes[indx]);
    r->setBirthTime(currentTime);
    r->setIsTip(true);
    r->setIsExtant(true);
    r->setIsExtinct(false);
    
    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);
    
    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(r);
    extantNodes.push_back(l);
    nodes.push_back(r);
    nodes.push_back(l);
    numExtant = (int)extantNodes.size();
    r->setIndx(numExtant - 2);
    l->setIndx(numExtant - 1);
    
}

void SpeciesTree::setBranchLengths(){
    double bl;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        bl = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(bl);
    }
}

void SpeciesTree::setPresentTime(double currentT){
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setDeathTime(currentT);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}

void SpeciesTree::setTreeTipNames(){
    unsigned extinctIt = 0;
    unsigned tipIt = 0;
    recTipNamer(getRoot(), extinctIt, tipIt);
}


void SpeciesTree::recTipNamer(Node *p, unsigned &extinctIndx, unsigned &tipIndx){
    if(p != NULL){
        std::stringstream tn;
        if(p->getIsTip()){
            if(p->getIsExtinct()){
                tn << p->getIndex();
                std::string name = "X" + tn.str();
                p->setName(name);
                
            }
            else{
                tn << p->getIndex();
                std::string name = "T" + tn.str();
                p->setName(name);
            }
        }
        else{
            recTipNamer(p->getLdes(), extinctIndx, tipIndx);
            recTipNamer(p->getRdes(), extinctIndx, tipIndx);
            
        }
    }
}

void SpeciesTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL)
            ss <<  p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ")";        }
    }
}



std::string SpeciesTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string spTree = ss.str();
    return spTree;
}

void SpeciesTree::setGSATipTreeFlags(){
    zeroAllFlags();
    numTotalTips = 0;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if((*it)->getIsTip()){
            numTotalTips++;
            (*it)->setFlag(1);
            
        }
        else{
            (*it)->setFlag(2);
        }
    }
    setSampleFromFlags();
}


//void SpeciesTree::setSampleFromFlags(){
//    int flag;
//    Node *q = NULL;
//    for(std::vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); p++){
//        if((*p)->getIsTip()){
//            flag = (*p)->getFlag();
//            q = (*p);
//            if(flag == 1){
//                while(q->getIsRoot() == false && flag < 2){
//                    q = q->getAnc();
//                    flag = q->getFlag();
//                    flag++;
//                    q->setFlag(flag);
//                    
//                }
//            }
//        }
//    }
//}

void SpeciesTree::popNodes(){
    nodes.clear();
    extantNodes.clear();
    recPopNodes(getRoot());
    
    int indx;
    for(std::vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); p++){
        indx = (int) (p - nodes.begin());
        (*p)->setIndx(indx);
    }
}

void SpeciesTree::recPopNodes(Node *p){
    if(p->getIsTip()){
        if(p->getIsExtant()){
            extantNodes.push_back(p);
            nodes.push_back(p);
        }
        else{
            nodes.push_back(p);
        }
    }
    else{
        nodes.push_back(p);
        recPopNodes(p->getLdes());
        recPopNodes(p->getRdes());
    }
}

void SpeciesTree::reconstructTreeFromGSASim(Node *oRoot){
    Node *n = new Node();
    unsigned tipCounter = extantStop;
    unsigned intNodeCounter = 0;
    reconstructLineageFromGSASim(n, oRoot, tipCounter, intNodeCounter);
    delete n;
}

void SpeciesTree::reconstructLineageFromGSASim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter){
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
                reconstructLineageFromGSASim(s1, prevN->getLdes(), tipCounter, intNodeCounter);
            if(prevN->getRdes()->getFlag() > 0)
                reconstructLineageFromGSASim(s1, prevN->getRdes(), tipCounter, intNodeCounter);
            
            
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
                reconstructLineageFromGSASim(currN, prevN->getLdes(), tipCounter, intNodeCounter);
            else
                reconstructLineageFromGSASim(currN, prevN->getRdes(), tipCounter, intNodeCounter);
        }
    }
}

std::map<int,double> SpeciesTree::getBirthTimesFromNodes(){
    int indx;
    double birthTime;
    std::map<int,double> birthTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        indx = (*it)->getIndex();
        birthTime = (*it)->getBirthTime();
        birthTimeMap.insert(std::pair<int,double>(indx, birthTime));
    }
    return birthTimeMap;
}

std::map<int,double> SpeciesTree::getDeathTimesFromNodes(){
    int indx;
    double deathTime;
    std::map<int,double> deathTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if(!((*it)->getIsExtant())){
            indx = (*it)->getIndex();
            deathTime = (*it)->getDeathTime();
            deathTimeMap.insert(std::pair<int,double>(indx, deathTime));
        }
    }
    return deathTimeMap;
}

std::pair<int,int> SpeciesTree::preorderTraversalStep(int indx){
    std::pair<int,int> sibs;
    sibs.first = nodes[indx]->getLdes()->getIndex();
    sibs.second = nodes[indx]->getRdes()->getIndex();
    return sibs;
}

int SpeciesTree::postOrderTraversalStep(int index){
    int d;
    d = nodes[index]->getAnc()->getIndex();
    return d;
}

bool SpeciesTree::macroEvent(int indx){
    bool isSpec;
    Node* n = nodes[indx];
    if(n->getIsTip())
        isSpec = false;
    else
        isSpec = true;
    return isSpec;
}
