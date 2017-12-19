 //
//  LocusTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/13/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "LocusTree.h"
#include <iostream>

LocusTree::LocusTree(MbRandom *p, unsigned nt, double stop, double gbr, double gdr, double lgtrate) : Tree(p, nt, 0.0){
    rando = p;
    numTaxa = 1;
    stopTime = stop;
    currentTime = 0.0;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtrate;
}

LocusTree::~LocusTree(){
    
}



void LocusTree::setNewLineageInfo(int indx, Node *r, Node *l){
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
    r->setIndx(extantNodes[indx]->getIndex());

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);
    l->setIndx(extantNodes[indx]->getIndex());

    
    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(r);
    extantNodes.push_back(l);
    nodes.push_back(r);
    nodes.push_back(l);
    numExtant = (int)extantNodes.size();
}

void LocusTree::lineageBirthEvent(int indx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    setNewLineageInfo(indx, right, sis);
}

void LocusTree::lineageDeathEvent(int indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void LocusTree::lineageTransferEvent(int indx){
    // TODO: write this puppy
    //first a birth event
    Node *donor, *rec;
    donor = new Node();
    rec = new Node();
    
    // donor keeps all the attributes  of the Node at extantNodes[indx]
    donor->setAnc(extantNodes[indx]);
    donor->setBirthTime(currentTime);
    donor->setIndx(extantNodes[indx]->getIndex());
    donor->setIsExtant(true);
    donor->setIsTip(true);
    donor->setIsExtinct(false);
    donor->setSib(extantNodes[indx]->getSib());
    donor->setLdes(NULL);
    donor->setRdes(NULL);
    
    
    //extantNodes[indx] goes extinct
    extantNodes[indx]->setLdes(rec);
    extantNodes[indx]->setRdes(donor);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setFlag(1);
    extantNodes[indx]->setIsExtant(false);
    

    // actual transfer event
    
    
    std::map<int,int> speciesIndx;
    unsigned randomSpeciesID;
    std::pair<int,int> recIndx;
    // need to draw a new ->getIndex
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        if((*it)->getIndex() != extantNodes[indx]->getIndex()){
            speciesIndx.insert(std::pair<int,int> (it - extantNodes.begin(), (*it)->getIndex()));
        }
    }

    randomSpeciesID = rando->uniformRv(0, (unsigned) speciesIndx.size() - 1);
    std::map<int,int>::iterator item = speciesIndx.begin();
    std::advance( item, randomSpeciesID );
    recIndx = *item;
    
    rec->setIndx(recIndx.second);
    rec->setBirthTime(currentTime);
    rec->setIsExtant(true);
    rec->setIsTip(true);
    rec->setIsExtinct(false);
    rec->setLdes(NULL);
    rec->setRdes(NULL);
    rec->setAnc(extantNodes[indx]);
    rec->setSib(extantNodes[recIndx.first]);
    
    speciesIndx.clear();
    
    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(rec);
    extantNodes.push_back(donor);
    nodes.push_back(rec);
    nodes.push_back(donor);
    numExtant = (int) extantNodes.size();
}

double LocusTree::getTimeToNextEvent(){
    double sumrt = geneBirthRate + geneDeathRate + transferRate;
    double returnTime = 0.0;
    returnTime = -log(rando->uniformRv()) / (stopTime * sumrt);
    currentTime += returnTime;
    return returnTime;
}

void LocusTree::ermEvent(double ct){
    double relBr = geneBirthRate / (geneDeathRate + geneBirthRate + transferRate);
    double relLGTr = transferRate / (geneBirthRate + geneDeathRate + transferRate) + relBr;
    double whichEvent = rando->uniformRv();
    unsigned long extantSize = extantNodes.size();
    unsigned nodeInd = rando->uniformRv(0, extantSize - 1);
    if(whichEvent < relBr){
        lineageBirthEvent(nodeInd);
    }
    else{
        if(whichEvent < relLGTr){
            if(numTaxa > 1 && extantSize > 1)
                lineageTransferEvent(nodeInd);
        }
        else
            lineageDeathEvent(nodeInd);
    }
}



int LocusTree::speciationEvent(int indx, double time){
    // indx is the index of the species that is to speciate at the input time
    Node *r, *l;
    int lociExtNodesIndx;
    int count = 0;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end();){
        lociExtNodesIndx = (*it)->getIndex();
        if(lociExtNodesIndx == indx){
            r = new Node();
            l = new Node();
 
            r->setLdes(NULL);
            r->setRdes(NULL);
            r->setSib(l);
            r->setAnc((*it));
            r->setBirthTime(time);
            r->setIsTip(true);
            r->setIsExtant(true);
            r->setIsExtinct(false);
            r->setIndx(-1);

            l->setLdes(NULL);
            l->setRdes(NULL);
            l->setSib(r);
            l->setAnc((*it));
            l->setBirthTime(time);
            l->setIsTip(true);
            l->setIsExtinct(false);
            l->setIsExtant(true);
            l->setIndx(-1);
            
            (*it)->setLdes(l);
            (*it)->setRdes(r);
            (*it)->setDeathTime(time);
            (*it)->setIsTip(false);
            (*it)->setIsExtant(false);
            nodes.push_back(r);
            nodes.push_back(l);
            it = extantNodes.erase(it);

            it = extantNodes.insert(it, r);
            it = extantNodes.insert(it, l);
            count += 2;
            numExtant = (int)extantNodes.size();
        }
        else
            ++it;
            
    }
    numTaxa++;
    return count;
}

void LocusTree::extinctionEvent(int indx, double time){
    // indx is the index of the species that is to go extinct at the input time
    int lociExtNodesIndx;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end();){
        lociExtNodesIndx = (*it)->getIndex();
        if(lociExtNodesIndx == indx){
            (*it)->setDeathTime(time);
            (*it)->setIsExtant(false);
            (*it)->setIsTip(true);
            (*it)->setIsExtinct(true);
            it = extantNodes.erase(it);
            numExtinct += 1;
            numExtant = (int) extantNodes.size();
        }
        else{
            it++;
        }
    }
    numTaxa--;
}

void LocusTree::setNewIndices(int indx, std::pair<int,int> sibs, int count){
    int lociExtNodesIndx;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end();){
        lociExtNodesIndx = (*it)->getIndex();
        if(lociExtNodesIndx == -1){
            (*it)->setIndx(sibs.first);
            (*it)->getSib()->setIndx(sibs.second);
            it += 2;
            count -= 2;
            if(count == 0)
                return;
        }
        else
            it++;
        
    }
}


void LocusTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL)
            ss << p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getLdes(), ss);
            ss << ":" << p->getLdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getRdes(), ss);
            ss << ":" << p->getRdes()->getBranchLength();
            ss << ")";
            
        }
    }
}

void LocusTree::setPresentTime(double currentT){
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setDeathTime(currentT);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}



std::string LocusTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string spTree = ss.str();
    return spTree;
}

void LocusTree::setTreeTipNames(){
    unsigned copyNumber = 0;
    bool recName = false;
    if(recName)
        recTipNamer(getRoot(), copyNumber);
    else{
        std::vector<int> copyNumberCounts;
        std::stringstream tn;
        std::string name;
        for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
            if((*it)->getIsTip()){
                if((*it)->getIsExtinct()){
                    tn << (*it)->getIndex();
                    copyNumberCounts.push_back((*it)->getIndex());
                    name = "X" + tn.str();
                    tn.clear();
                    tn.str(std::string());
                    tn << std::count(copyNumberCounts.begin(),copyNumberCounts.end(), (*it)->getIndex());
                    name += "_" + tn.str();
                    (*it)->setName(name);
                    tn.clear();
                    tn.str(std::string());


                }
                else{
                    tn << (*it)->getIndex();
                    copyNumberCounts.push_back((*it)->getIndex());
                    name = "T" + tn.str();
                    tn.clear();
                    tn.str(std::string());
                    tn << std::count(copyNumberCounts.begin(),copyNumberCounts.end(),  (*it)->getIndex());
                    name += "_" + tn.str();
                    (*it)->setName(name);
                    tn.clear();
                    tn.str(std::string());

                }
            }
            else{
                if((*it)->getFlag() == 1){
                    tn << (*it)->getIndex();
                    copyNumberCounts.push_back((*it)->getIndex());
                    name = "TR" + tn.str();
                    tn.clear();
                    tn.str(std::string());
                    tn << std::count(copyNumberCounts.begin(),copyNumberCounts.end(),  (*it)->getIndex());
                    name += "_" + tn.str();
                    (*it)->setName(name);
                    tn.clear();
                    tn.str(std::string());

                }
            
            }
        }
    }
}

// NOTE: this names tips but doesn't have unique tip names
void LocusTree::recTipNamer(Node *p, unsigned &copyNumber){
    if(p != NULL){
        std::stringstream tn;
        if(p->getIsTip()){
            if(p->getIsExtinct()){
                tn << p->getIndex();
                std::string name = "X" + tn.str();
                tn << copyNumber;
                name += "_" + tn.str();
                p->setName(name);
                copyNumber++;
                
            }
            else{
                tn << p->getIndex();
                std::string name = "T" + tn.str();
                tn << copyNumber;
                name += "_" + tn.str();
                p->setName(name);
                copyNumber++;
            }
        }
        else{
            recTipNamer(p->getLdes(), copyNumber);
            copyNumber = 0;
            recTipNamer(p->getRdes(), copyNumber);
            
        }
    }
}


void LocusTree::setBranchLengths(){
    double brlen;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        brlen = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(brlen);  
    }
}