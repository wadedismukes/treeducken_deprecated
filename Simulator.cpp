//
//  Simulator.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Simulator.h"

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

Simulator::~Simulator(){
    if(spTree != 0){
        delete spTree;
    }
    if(gsaTrees.size() != 0)
        gsaTrees.clear();
    if(lociTrees.size() != 0)
       lociTrees.clear();
}

void Simulator::initializeSim(){
    // only simulates speciestrees, other models to come
    if(simType == 1){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
    }
    else if(simType == 2){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
        lociTrees.resize(numLoci, new LocusTree(rando, numTaxaToSim, currentSimTime, geneBirthRate, geneDeathRate, transferRate)); // the 0.0 is the stopping time which will be set after the gsa CBDP sim of the species tree
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
    double stopTime = currentSimTime;
    double eventTime;
    bool isSpeciation;
    int countCopies;
    
    std::map<int,double> speciesBirthTimes = spTree->getBirthTimesFromNodes();
    std::map<int,double> speciesDeathTimes = spTree->getDeathTimesFromNodes();
    std::set<int> contempSpecies;
    std::pair<int, int> sibs;
    
    Node* spRoot = spTree->getRoot();
    for(int i = 0; i < numLoci; i++){
        lociTrees[i]->setStopTime(stopTime);
        currentSimTime = 0;
        contempSpecies.clear();
        contempSpecies.insert(spRoot->getIndex());
        while(currentSimTime < stopTime){
            eventTime = lociTrees[i]->getTimeToNextEvent();
            currentSimTime += eventTime;
            for(std::set<int>::iterator it = contempSpecies.begin(); it != contempSpecies.end();){
                if(currentSimTime > speciesDeathTimes[(*it)]){
               //     std::cout << speciesDeathTimes[(*it)] << std::endl;
                    isSpeciation = spTree->macroEvent(*it);
                    if(isSpeciation){
                        countCopies = lociTrees[i]->speciationEvent((*it), speciesDeathTimes[(*it)]);
                        sibs = spTree->preorderTraversalStep(*it);
                        lociTrees[i]->setNewIndices(*it, sibs, countCopies);
                        it = contempSpecies.erase(it);
                        it = contempSpecies.insert( it, sibs.first);
                        it = contempSpecies.insert( it, sibs.second);
                        --it;
                    }
                    else{
                        if(!(spTree->getIsExtantFromIndx(*it))){
                            lociTrees[i]->extinctionEvent(*it, speciesDeathTimes[(*it)]);
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
                
                if(lociTrees[i]->getNumExtant() < 1){
                    treesComplete = false;
                    return treesComplete;
                }
            }
            if(currentSimTime >= stopTime)
                currentSimTime = stopTime;
            else
                lociTrees[i]->ermEvent(currentSimTime);
            if(lociTrees[i]->getNumExtant() < 1){
                treesComplete = false;
                return treesComplete;
            }



        }

        lociTrees[i]->setPresentTime(stopTime);
    }
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
    return lociTrees[0]->printNewickTree();
}

