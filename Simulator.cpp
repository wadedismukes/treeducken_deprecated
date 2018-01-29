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

Simulator::Simulator(MbRandom *p, unsigned ntax, double lambda, double mu, double rho, unsigned numLociToSim, double gbr, double gdr, double lgtr, unsigned ipp, unsigned Ne, double genTime)
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
    generationTime = genTime;
    
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
    }
    else if(simType == 3){
        spTree = new SpeciesTree(rando, numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
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
    lociTree = new LocusTree(rando, numTaxaToSim, currentSimTime, geneBirthRate, geneDeathRate, transferRate);

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
        for(std::set<int>::iterator it = contempSpecies.begin(); it != contempSpecies.end();){
            if(currentSimTime > speciesDeathTimes[(*it)]){
                isSpeciation = spTree->macroEvent((*it));
                if(isSpeciation){
              //      std::cout << "species death time of species " << (*it) << " is " << speciesDeathTimes[(*it)] << std::endl;
                    sibs = spTree->preorderTraversalStep(*it);
                    countCopies = lociTree->speciationEvent((*it), speciesDeathTimes[(*it)], sibs);
            //        std::cout << sibs.first << "," << sibs.second << std::endl;
                    it = contempSpecies.erase(it);
                    it = contempSpecies.insert( it, sibs.second);
                    ++it;
                    it = contempSpecies.insert( it, sibs.first);
                }
                else{
                    if(!(spTree->getIsExtantFromIndx(*it))){
                    //    std::cout << "extinction of species " << (*it)   << " at " << speciesDeathTimes[(*it)] << std::endl;
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
          //  std::cout << "current sim time before event: " << currentSimTime << std::endl;
            lociTree->ermEvent(currentSimTime);
        }



    }
    lociTree->setPresentTime(currentSimTime);
    treesComplete = true;

    return treesComplete;
}

bool Simulator::simSpeciesLociTrees(){
    bool good = false;
    bool spGood = false;
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


std::set<double, std::greater<double>> Simulator::getEpochs(std::multimap<int,double> birthMap){
    
    std::set<double, std::greater<double>> epochs;
    
    for(std::multimap<int, double>::iterator it = birthMap.begin(); it != birthMap.end(); ++it){
        epochs.insert((*it).second);
    }
    
    
    return epochs;
}


bool Simulator::coalescentSim(){
    bool treeGood = false;
    bool reachedEnd = true;
    int ancIndx;
    geneTree = new GeneTree(rando, numTaxaToSim, indPerPop, popSize, generationTime);

    currentSimTime *= (generationTime);
    std::map<int,int> spToLo;
    
    // map with keys as indices of lociTree nodes vector, values birth times
    std::multimap<int, double> locusBirthMap = lociTree->getBirthTimesFromNodes();
    // scale to generations
    locusBirthMap = geneTree->rescaleTimes(locusBirthMap);
    
    // map with keys as indices of lociTree nodes vector, values death times
    std::multimap<int, double> locusDeathMap = lociTree->getDeathTimesFromNodes();
    // scale to generations
    locusDeathMap = geneTree->rescaleTimes(locusDeathMap);
    
    // map with keys as indices of lociTree nodes vector, values death times, this includes ONLY extinct loci
    std::multimap<int,double> deadSpeciesStartTimes = lociTree->getDeathTimesFromExtinctNodes();
    // scale to generations
    deadSpeciesStartTimes = geneTree->rescaleTimes(deadSpeciesStartTimes);
    
    // unordered set containing loci that are currently alive
    std::unordered_set<int> contempLoci;
    // get all time slices to simulate through
    std::set<double, std::greater<double>> epochs = getEpochs(locusDeathMap);
    
    currentSimTime = *(epochs.begin());
    epochs.erase(epochs.begin());
    epochs.insert(epochs.end(), 0.0);
    // need to populate extantNodes with extant loci (indPerPop for each loci)
    contempLoci = lociTree->getExtantLoci();
    // set Node values in geneTree
    geneTree->initializeTree(contempLoci, locusDeathMap);
    std::unordered_set<int>::iterator contempLociEnd;
    
    
    for(std::set<double, std::greater<double>>::iterator it = epochs.begin(); it != epochs.end(); ++it){
        
        contempLociEnd = contempLoci.end();
        if(*it != 0.0){
            for(std::unordered_set<int>::iterator locIt = contempLoci.begin(); locIt != contempLociEnd; ){
                                 
                ancIndx = lociTree->postOrderTraversalStep(*locIt);
                reachedEnd = geneTree->censorCoalescentProcess(currentSimTime, *it, *locIt, ancIndx);
                
                if(!(reachedEnd)){
                    locIt = contempLoci.erase(locIt);
                    contempLoci.insert(locIt, ancIndx);
                }
                else{
                    ++locIt;
                }
            }
            
            for(std::multimap<int,double>::iterator temp = deadSpeciesStartTimes.begin(); temp != deadSpeciesStartTimes.end(); ++temp){
                if(temp->second == *it){
                    std::cout << "epoch time: " << *it << std::endl;
                    std::cout << "added locus " << temp->first << std::endl;
                    geneTree->addExtinctSpecies(temp->second, temp->first);
                    contempLoci.insert(temp->first);
                }
            }
            currentSimTime = *it;
        }
        else{
            currentSimTime = 0.0;
            geneTree->rootCoalescentProcess(currentSimTime);
            treeGood = true;
            spToLo = lociTree->getLocusToSpeciesMap();
            geneTree->setIndicesBySpecies(spToLo);
            break;
        }

    }

    return treeGood;
}

bool Simulator::simThreeTree(){
    bool gGood = false;
    bool spGood = false;
    bool loGood = false;
    
    while(!gGood){
        while(!loGood){
            while(!spGood){
                spGood = gsaBDSim();
            }
            loGood = bdsaBDSim();
        }
        gGood = coalescentSim();
    }
    
    return gGood;
}


std::string Simulator::printGeneTreeNewick(){
    std::string newickTree;
    newickTree = geneTree->printNewickTree();
    return newickTree;
}