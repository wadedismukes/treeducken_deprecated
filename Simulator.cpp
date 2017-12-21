//
//  Simulator.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Simulator.h"
#include <iostream>

Simulator::Simulator(MbRandom *p, unsigned nt, double lambda, double mu, double rho)
{
    simType = 1;
    currentSimTime = 0.0;
    rando = p;
    numTaxaToSim = nt;
    gsaStop = 100*nt;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    
    
    numLoci = 0;
    numGenes = 0;
    geneBirthRate = 0.0;
    geneDeathRate = 0.0;
    transferRate = 0.0;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0;
    
}


Simulator::Simulator(MbRandom *p, unsigned ntax, double lambda, double mu, double rho, unsigned numLociToSim, double gbr, double gdr, double lgtr)
{
    simType = 2;
    currentSimTime = 0.0;
    rando = p;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    
    numLoci = numLociToSim;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0;
    
    
}

Simulator::Simulator(MbRandom *p, unsigned ntax, double lambda, double mu, double rho, unsigned numLociToSim, double gbr, double gdr, double lgtr, unsigned ipp, unsigned Ne)
{
    simType = 3;
    currentSimTime = 0.0;
    rando = p;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    
    numLoci = numLociToSim;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = ipp;
    popSize = Ne;
    
}

Simulator::~Simulator(){
    if(spTree != 0){
        delete spTree;
    }
    if(gsaTrees.size() != 0)
        gsaTrees.clear();
    if(lociTree != 0)
        delete lociTree;
    if(geneTree != 0)
        delete geneTree;
        
}

void Simulator::initializeSim(){
    // only simulates speciestrees, other models to come
    if(simType == 1){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
    }
    else if(simType == 2){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
        lociTree = new LocusTree(rando, numTaxaToSim, currentSimTime, geneBirthRate, geneDeathRate, transferRate);
    }
    else if(simType == 3){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
        lociTree = new LocusTree(rando, numTaxaToSim, currentSimTime, geneBirthRate, geneDeathRate, transferRate);
        geneTree = new GeneTree(rando, numTaxaToSim, indPerPop, popSize);

        // geneTrees.resize(numLoci, new GeneTree(rando, numTaxaToSim, currentSimTime, indPerPop, popSize);
    }
}


/*
 Below is the machinery to use GSA sampling (Hartmann 2010) to simulate a species tree.
 Much of this code is modified from FossilGen (written by Tracy Heath)
 */
bool Simulator::gsaBDSim(){
    double timeIntv, sampTime;
    bool treeComplete = false;
    this->initializeSim();
    double eventTime;
    
    while(gsaCheckStop()){
        eventTime = spTree->getTimeToNextEvent();
        currentSimTime += eventTime;
        spTree->ermEvent(currentSimTime);
        if(spTree->getNumExtant() < 1){
            treeComplete = false;
            return treeComplete;
        }
        else if(spTree->getNumExtant() == numTaxaToSim){
            timeIntv = spTree->getTimeToNextEvent();
            sampTime = rando->uniformRv(0, timeIntv) + currentSimTime;
            spTree->setPresentTime(sampTime);
            processGSASim();
        }
        
    }
    unsigned gsaRandomTreeID = rando->uniformRv(0, (unsigned) gsaTrees.size() - 1);
    spTree = gsaTrees[gsaRandomTreeID];
    processSpTreeSim();
    spTree->setBranchLengths();
    spTree->setTreeTipNames();
    currentSimTime = spTree->getCurrentTimeFromExtant();
    //spTree->setPresentTime(currentSimTime);
    
    treeComplete = true;
    
    return treeComplete;
}


bool Simulator::gsaCheckStop(){
  
  bool keepSimulating = true;
  
  if(spTree->getNumExtant() >= gsaStop){
      keepSimulating = false;
  }
  
  return keepSimulating;
}

void Simulator::processGSASim(){
    SpeciesTree *tt = new SpeciesTree(rando, numTaxaToSim + spTree->getNumExtinct());
    this->prepGSATreeForReconstruction();
    Node *simRoot = spTree->getRoot();
    tt->reconstructTreeFromGSASim(simRoot);
    gsaTrees.push_back(tt);
    
}

void Simulator::processSpTreeSim(){
    spTree->setSpeciationRate(speciationRate);
    spTree->setExtinctionRate(extinctionRate);
    spTree->popNodes();
}

void Simulator::prepGSATreeForReconstruction(){
    spTree->setGSATipTreeFlags();
}


bool Simulator::simSpeciesTree(){
    bool good = false;
    while(!good){
        good = gsaBDSim();
    }
    return good;
}

std::string Simulator::printSpeciesTreeNewick(){
    return spTree->printNewickTree();
}

bool Simulator::bdsaBDSim(){
    bool treesComplete = false;
    double stopTime = spTree->getCurrentTimeFromExtant();
    double eventTime;
    bool isSpeciation;
    int countCopies;
    
    std::map<int,double> speciesBirthTimes = spTree->getBirthTimesFromNodes();
    std::map<int,double> speciesDeathTimes = spTree->getDeathTimesFromNodes();
    std::set<int> contempSpecies;
    std::pair<int, int> sibs;
    
    Node* spRoot = spTree->getRoot();
    lociTree->setStopTime(stopTime);
    currentSimTime = 0;
    
    if(!(contempSpecies.empty()))
        contempSpecies.clear();
    contempSpecies.insert(spRoot->getIndex());
    
    while(currentSimTime < stopTime){
        eventTime = lociTree->getTimeToNextEvent();
        currentSimTime += eventTime;
        std::cout << "simTime" << currentSimTime << std::endl;
        for(std::set<int>::iterator it = contempSpecies.begin(); it != contempSpecies.end();){
            if(currentSimTime > speciesDeathTimes[(*it)]){
                // std::cout << speciesDeathTimes[(*it)] << std::endl;
                isSpeciation = spTree->macroEvent((*it));
                if(isSpeciation){
                    std::cout << "species death time of species " << (*it) << " is " << speciesDeathTimes[(*it)] << std::endl;
                    sibs = spTree->preorderTraversalStep(*it);
                    countCopies = lociTree->speciationEvent((*it), speciesDeathTimes[(*it)], sibs);
                    std::cout << sibs.first << "," << sibs.second << std::endl;
  //                  lociTrees[i]->setNewIndices(*it, sibs, countCopies);
                    it = contempSpecies.erase(it);
                    it = contempSpecies.insert( it, sibs.second);
                    it++;
                    it = contempSpecies.insert( it, sibs.first);
                }
                else{
                    if(!(spTree->getIsExtantFromIndx(*it))){
                        std::cout << "extinction of species " << (*it)   << " at " << speciesDeathTimes[(*it)] << std::endl;
                        lociTree->extinctionEvent(*it, speciesDeathTimes[(*it)]);
                        it = contempSpecies.erase(it);
                    }
                    else{
                        ++it;
                    }
                }
            }
            else{
                ++it;
            }
            
            if(lociTree->getNumExtant() < 1){
                treesComplete = false;
                return treesComplete;
            }
            
        }
        
        if(currentSimTime >= stopTime){
            currentSimTime = stopTime;
            lociTree->setCurrentTime(stopTime);
        }
        else{
            std::cout << "current sim time before event: " << currentSimTime << std::endl;
            lociTree->ermEvent(currentSimTime);
        }



    }
    std::cout << "currentSimTime at end: " << currentSimTime << std::endl;
    std::cout << "currentTime at end: " << lociTree->getCurrentTime() << std::endl;
    lociTree->setPresentTime(currentSimTime);
    treesComplete = true;

    return treesComplete;
}

bool Simulator::simSpeciesLociTrees(){
    bool good, spGood = false;
    while(!good){
        while(!spGood){
            spGood = gsaBDSim();
        }
        good = bdsaBDSim();
    }
    
    return good;
}

std::string Simulator::printLocusTreeNewick(){
    std::string newickTree;
    newickTree = lociTree->printNewickTree();
    return newickTree;
}

