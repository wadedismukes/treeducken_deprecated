#include "Simulator.h"
#include <iostream>

/**
 * Constructor of Simulator class for full three-tree model
 * @param p RNG seed from MrBayes code
 * @param ntax Number of taxa to simulate
 * @param lambda Rate of birth as double
 * @param mu Rate of death as double
 * @param rho Sampling Fraction (currently not implemented)
 * @param numLociToSim Number of loci to simulate per species tre
 * @param gbr Gene birth rate as double
 * @param gdr Gene death rate as double
 * @param lgtr Lateral gene transfer rate
 * @param ipp Individuals per locus tree lineage to perform coalescent process on
 * @param Ne Total population size of each locus tree lineage
 * @param genTime Scaling factor for converting from coalescent units to absolute time entered as generations per unit time
 * @param ng Number of gene trees to simulate for each locus tree lineage
 * @param og Outgroup fraction, pastes a fake outgroup on tree with branch length set to some fraction of total tree depth
 * @param ts Scaling factor for tree entered as a double to convert trees to some length (e.g. 1.0)
 * @param sout Bool for turning off stdout (speeds up the program)
 */
Simulator::Simulator(MbRandom *p,
                     unsigned ntax,
                     double lambda,
                     double mu,
                     double rho,
                     unsigned numLociToSim,
                     double gbr,
                     double gdr,
                     double lgtr,
                     unsigned ipp,
                     unsigned Ne,
                     double genTime,
                     int ng,
                     double og,
                     double ts,
                     bool sout)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 3;
    currentSimTime = 0.0;
    rando = p;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    numLoci = numLociToSim;
    numGenes = ng;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = ipp;
    popSize = Ne;
    printSOUT = sout;
    generationTime = genTime;
    outgroupFrac = og;
    geneTrees.resize(numLoci);
    treeScale = ts;
    propDuplicate = -1;
}
/**
 * Destructor for Simulator classes
 */
Simulator::~Simulator(){
    for(auto & gsaTree : gsaTrees){
        delete gsaTree;
    }
    gsaTrees.clear();
    int i = 0;
    for(auto & locusTree : locusTrees){
        delete locusTree;
        for(auto & q : geneTrees[i]){
            delete q;
        }
        geneTrees[i].clear();
        ++i;
    }
    locusTrees.clear();
    
        
}


/**
 * Generalized Sampling Algorithm for generating birth-death trees of the correct length
 *
 * @details Below is the machinery to use GSA sampling (Hartmann 2010) to simulate a species tree. Much of this code is modified from FossilGen (written by Tracy Heath)
 */
bool Simulator::gsaBDSim(){
    double timeInterval, sampTime;
    bool treeComplete;
    auto st = SpeciesTree(rando, numTaxaToSim, speciationRate, extinctionRate);
    spTree = &st;
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
            timeInterval = spTree->getTimeToNextEvent();
            sampTime = rando->uniformRv(0, timeInterval) + currentSimTime;
            spTree->setPresentTime(sampTime);
            processGSASim();
        }
        
    }
    unsigned gsaRandomTreeID = rando->uniformRv(0, (unsigned) gsaTrees.size() - 1);
    // delete spTree;
    spTree = gsaTrees[gsaRandomTreeID];
    processSpTreeSim();
    spTree->setBranchLengths();
    spTree->setTreeTipNames();
    currentSimTime = spTree->getCurrentTimeFromExtant();
    if(treeScale > 0.0){
        spTree->scaleTree(treeScale, currentSimTime);
        currentSimTime = treeScale;
    }

    treeComplete = true;
    
    return treeComplete;
}

/**
 * Function for checking whether the simulation has reached the stopping point
 * @return Bool for whether to keepSimulating or not
 */

bool Simulator::gsaCheckStop(){

    bool keepSimulating;
    keepSimulating = true;
  
  if(spTree->getNumExtant() >= gsaStop){
      keepSimulating = false;
  }
  
  return keepSimulating;
}

/**
 * Function for completing the simulation
 * @details this prunes the desired species tree from the larger simulated tree according to Hartmann et al. 2010
 */
void Simulator::processGSASim(){
    auto *tt = new SpeciesTree(rando, numTaxaToSim + spTree->getNumExtinct());
    this->prepGSATreeForReconstruction();
    Node *simRoot = spTree->getRoot();
    tt->setRoot(simRoot);
    tt->reconstructTreeFromGSASim(simRoot);
    gsaTrees.push_back(tt);    
}

/**
 * Function for processing the final result of a species tree simulation
 * @details This is performed after the tree is rebuilt using the recursive function reconstructTreeFromGSASim, repopulates the nodes vector in SpeciesTree class
 *
 * @see reconstructTreeFromGSASim()
 */
void Simulator::processSpTreeSim(){
    spTree->setSpeciationRate(speciationRate);
    spTree->setExtinctionRate(extinctionRate);
    spTree->popNodes();

}

/**
 * Sets flags from the SpeciesTree class for pruning
 *
 * @see reconstructTreeFromGSASim()
 */
void Simulator::prepGSATreeForReconstruction(){
    spTree->setGSATipTreeFlags();
}

/**
 * Function for calling simulation function for a species tree under a Moran model
 *
 * @details Currently does not work
 *
 * @return A bool indicating whether the simulation completed correctly
 */
bool Simulator::simMoranSpeciesTree(){
    auto simulationComplete = false;
    while(!simulationComplete){
        simulationComplete = moranSpeciesSim();
    }
    return simulationComplete;
}

/**
 * Simulates a species tree under a Moran process
 *
 * @todo fix this function
 * @return a bool passed to simMoranSpeciesTree for ensuring a complete tree is simulated
 */
bool Simulator::moranSpeciesSim(){
    bool treeComplete;
    SpeciesTree st = SpeciesTree(rando, numTaxaToSim, speciationRate, extinctionRate);
    spTree = &st;
    spTree->initializeMoranProcess(numTaxaToSim);
    double eventTime;
    
    while(moranCheckStop()){
        eventTime = spTree->getTimeToNextEventMoran();
        currentSimTime += eventTime;
        spTree->moranEvent(currentSimTime);
        if(spTree->getNumExtant() < 1){
            treeComplete = false;
            return treeComplete;
        }
    }
    // processing will make a tree with extant and non extant tips
    //spTree->setBranchLengths();
    //spTree->setTreeTipNames();
    if(treeScale > 0.0){
        spTree->scaleTree(treeScale, currentSimTime);
        currentSimTime = treeScale;
    }

    treeComplete = true;
    
    return treeComplete;
}

/**
 * Checks if simulation is complete
 * @return True to keep simulating, false otherwise
 */
bool Simulator::moranCheckStop(){

    auto keepSimulating = true;
  
  if((spTree->getNumExtant() * 2) <= currentSimTime){
      keepSimulating = false;
  }
  
  return keepSimulating;
}

/**
 * Calls simulation function and ensures that function provides a good simulation
 * @return a bool if the simulation is completed succesfully true, false otherwise
 */
bool Simulator::simSpeciesTree(){
    bool good = false;
    while(!good){
        good = gsaBDSim();
    }
    if(outgroupFrac > 0.0)
        this->graftOutgroup(spTree, spTree->getTreeDepth());
    return good;
}

/**
 * Prints a Newick tree with only extant species on it pruning the remaining taxa
 * @return A Newick string
 */
std::string Simulator::printExtSpeciesTreeNewick(){
    auto *tt = new SpeciesTree(rando, numTaxaToSim);
    spTree->getRootFromFlags(false);
    if(outgroupFrac > 0.0){
        tt->setOutgroup(spTree->getOutgroup());
        tt->setRoot(spTree->getOutgroup()->getAnc());
    }
    else{
        tt->setRoot(spTree->getExtantRoot());
    }
    tt->setExtantRoot(tt->getRoot());
    tt->reconstructTreeFromSim(spTree->getRoot());
    std::string newickTree = tt->printExtNewickTree();
    delete tt;
    return newickTree;
}

/**
 * Wrapper for printNewickTree of Tree class
 * @return A Newick string showing the structure of the SpeciesTree class held in spTree
 */
std::string Simulator::printSpeciesTreeNewick(){
    return spTree->printNewickTree();
}

/**
 * Simulates a locus tree storing the tree stucture in the LocusTree class
 * @return a bool indicating if all the simulated trees simulated to completion
 */
bool Simulator::bdsaBDSim(){
    bool treesComplete;
    double stopTime = spTree->getCurrentTimeFromExtant();
    double eventTime;
    bool isSpeciation;
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
        for(auto it = contempSpecies.begin(); it != contempSpecies.end();){
            if(currentSimTime > speciesDeathTimes[(*it)]){
                isSpeciation = spTree->macroEvent((*it));
                if(isSpeciation){
                    sibs = spTree->tipwiseStep(*it);
                    lociTree->speciationEvent((*it), speciesDeathTimes[(*it)], sibs);
                    it = contempSpecies.erase(it);
                    it = contempSpecies.insert( it, sibs.second);
                    ++it;
                    it = contempSpecies.insert( it, sibs.first);
                }
                else{
                    if(!(spTree->getIsExtantFromIndx(*it))){
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
            lociTree->ermEvent(currentSimTime);
        }



    }
    lociTree->setPresentTime(currentSimTime);
    treesComplete = true;

    return treesComplete;
}

/**
 * Function that calls gsaBDSim and then for each of loci bdsaSim to generate numLoci classes of LocusTree
 * @return a bool indicating if both completed
 */
bool Simulator::simSpeciesLociTrees(){
    bool good = false;
    bool spGood = false;
    for(int i = 0; i < numLoci; i++){
        while(!good){
            while(!spGood){
                spGood = gsaBDSim();
            }
            if(outgroupFrac > 0.0)
                this->graftOutgroup(spTree, spTree->getTreeDepth());
            if(printSOUT)
                std::cout << "Simulating loci #" <<  i + 1 << std::endl;
            good = bdsaBDSim();
        }
        if(outgroupFrac > 0.0)
            this->graftOutgroup(lociTree, lociTree->getTreeDepth());

        locusTrees.push_back(lociTree);

        good = false;
    }
    return good;
}

/**
 * Prints out a locus tree at index i from the vector of class LocusTree a
 * @param i index of the locus tree to be printed out
 * @return a Newick string of the locus tree at index i of the vector of class LocusTree
 */
std::string Simulator::printLocusTreeNewick(int i){
    std::string newickTree;
    auto it = locusTrees.begin();
    std::advance(it, i);
    newickTree = (*it)->printNewickTree();
    return newickTree;
}

/**
 * Function to create a set of sorted doubles of epochs
 * @details epochs are defined by branching points on the LocusTree and extinction events of nodes dying before present on the locus tree. sorted in reverse order
 * @return A set containing the epoch times of a LocusTree stored in lociTree
 */
std::set<double, std::greater<double> > Simulator::getEpochs(){
    std::set<double, std::greater<double> > epochs;
    std::vector<Node*> lociTreeNodes = lociTree->getNodes();
    std::vector<Node *, std::allocator<Node *>>::iterator it;
    for(it = lociTreeNodes.begin(); it != lociTreeNodes.end(); ++it){
        if(!((*it)->getIsExtinct())){
            if((*it)->getIsTip())
                epochs.insert((*it)->getDeathTime());
            epochs.insert((*it)->getBirthTime());
        }
        else
            epochs.insert((*it)->getDeathTime());
    }
    return epochs;
}

/**
 * Simulates a censored coalescent process and stores resulting coalescent tree in geneTree
 * @return A bool indicating if the censored coalescent process finished
 */
bool Simulator::coalescentSim(){
    bool treeGood = false;
    geneTree = new GeneTree(rando, numTaxaToSim, indPerPop, popSize, generationTime);

    std::map<int,int> spToLo;

    int ancIndx;
    int epochCount = 0;

    double stopTime, stopTimeEpoch, stopTimeLoci;
    bool deathCheck;
    bool allCoalesced;
    bool is_ext;

    std::set<double, std::greater<double> > epochs = getEpochs();
    int numEpochs = (int) epochs.size();
    std::set<int> extinctFolks = lociTree->getExtLociIndx();
    std::set<int> coalescentBounds = lociTree->getCoalBounds();
    std::vector< std::vector<int> > contempLoci = lociTree->getExtantLoci(epochs);
    std::map<int, double> stopTimes = lociTree->getBirthTimesFromNodes();
    geneTree->initializeTree(contempLoci, *(epochs.begin()));
    if(outgroupFrac != 0.0)
        contempLoci[0].pop_back();
    std::set<int>::iterator extFolksIt;

    std::set<double, std::greater<double>, std::allocator<double>>::iterator epIter;
    for(epIter = epochs.begin(); epIter != epochs.end(); ++epIter){
        currentSimTime = *epIter;
        if(epochCount != numEpochs - 1){
            epIter = std::next(epIter, 1);
            stopTimeEpoch = *epIter;
            for(int j = 0; j < contempLoci[epochCount].size(); ++j){
                extFolksIt = extinctFolks.find(contempLoci[epochCount][j]);
                is_ext = (extFolksIt != extinctFolks.end());
                if(is_ext){
                    geneTree->addExtinctSpecies(currentSimTime, contempLoci[epochCount][j]);
                    extinctFolks.erase(extFolksIt);
                }
                stopTimeLoci = stopTimes[contempLoci[epochCount][j]];
                
                if(stopTimeLoci > stopTimeEpoch){
                    stopTime = stopTimeLoci;
                    deathCheck = true;
                }
                else{
                    stopTime = stopTimeEpoch;
                    deathCheck = false;
                }

                ancIndx = lociTree->postOrderTraversalStep(contempLoci[epochCount][j]);
                allCoalesced = geneTree->censorCoalescentProcess(currentSimTime, stopTime, contempLoci[epochCount][j], ancIndx, deathCheck);
                
                
                // if all coalesced remove that loci from the matrix of loci
                if(allCoalesced){
                    int check = contempLoci[epochCount][j];
                    for(int k = epochCount + 1; k < numEpochs; k++){
                        for(int m = 0; m < contempLoci[k].size(); ++m){
                            if(contempLoci[k][m] == check){
                                contempLoci[k].erase(contempLoci[k].begin() + m);
                                break;
                            }
                        }
                    }
                }
            }
            epIter = std::prev(epIter, 1); 
        }
        else{
            // finish coalescing
            geneTree->rootCoalescentProcess(currentSimTime, outgroupFrac);
            treeGood = true;
        }
        epochCount++;
    }

    spToLo = lociTree->getLocusToSpeciesMap();
    geneTree->setIndicesBySpecies(spToLo);
    return treeGood;
}

/**
 * Wrapper for a full simulation of the three-tree model
 * @details First simulates a species tree, for each species tree numLoci LocusTree classes are simulated, for each LocusTree class numGene GeneTree classes are simulated
 * @return a bool indicating that the simulation completed
 */
bool Simulator::simThreeTree(){
    bool gGood = false;
    bool spGood = false;
    bool loGood = false;
    while(!spGood){
        spGood = gsaBDSim();

    }
    for(int i = 0; i < numLoci; i++){
        while(!loGood){
            if(printSOUT)
                std::cout << "Simulating loci # " <<  i + 1 << std::endl;
            loGood = bdsaBDSim();
        }
        if(outgroupFrac > 0.0){
            this->graftOutgroup(lociTree, lociTree->getTreeDepth());
        }
        for(int j = 0; j < numGenes; j++){
            while(!gGood){
                if(printSOUT)
                    std::cout << "Simulating gene # " <<  j + 1 << " of loci # " << i + 1 << std::endl;
                gGood = coalescentSim();
            }
            geneTrees[i].push_back(geneTree);

            gGood = false;
        }
        locusTrees.push_back(lociTree);
        loGood = false;
    }
    if(outgroupFrac > 0.0)
        this->graftOutgroup(spTree, spTree->getTreeDepth());

    return gGood;
}

/**
 * Prints Newick string for GeneTree j found in LocusTree i of simulator class
 * @param i index of LocusTree in Simulator class
 * @param j index of GeneTree in Simulator class
 * @return
 */
std::string Simulator::printGeneTreeNewick(int i, int j){
    std::string newickTree;
    newickTree = geneTrees[i][j]->printNewickTree();
    return newickTree;
}
/**
 * Prints Newick string containing only extant taxa for gene tree j found in LocusTree i of Simulator class
 * @param i
 * @param j
 * @return
 */
std::string Simulator::printExtantGeneTreeNewick(int i, int j){
    std::string newickTree;
    if(geneTrees[i][j]->getExtantNodes().size() > 1){
        auto *tt = new GeneTree(rando, numTaxaToSim, indPerPop, popSize, generationTime);
        geneTrees[i][j]->getRootFromFlags(true);
        if(outgroupFrac > 0.0){
            tt->setOutgroup(geneTrees[i][j]->getOutgroup());
            tt->setRoot(geneTrees[i][j]->getOutgroup()->getAnc());
        }
        else
            tt->setRoot(geneTrees[i][j]->getExtantRoot());
        tt->setExtantRoot(geneTrees[i][j]->getExtantRoot());
        tt->reconstructTreeFromSim(geneTrees[i][j]->getRoot());
        newickTree = tt->printExtantNewickTree();
        delete tt;
    }
    else{
        newickTree = ";";
    }
    return newickTree;
}

/**
 * Function graft an outgroup onto a tree
 *
 * @param tr tree to be grafted
 * @param trDepth Depth of the tree stored at *tr
 */
void Simulator::graftOutgroup(Tree *tr, double trDepth){
    auto *rootNode = new Node();
    Node *currentRoot = tr->getRoot();
    rootNode->setBirthTime(currentRoot->getBirthTime());
    auto *outgroupNode = new Node();
    tr->rescaleTreeByOutgroupFrac(outgroupFrac, trDepth);
    double tipTime = tr->getEndTime();
    tr->setNewRootInfo(rootNode, outgroupNode, currentRoot, tipTime);
}

/**
 * Wrapper for simulating locus and then gene trees
 *
 *
 * @return A bool indicating a completed simulation
 */
bool Simulator::simLocusGeneTrees(){
    bool loGood = false;
    bool gGood = false;
    for(int i = 0; i < numLoci; i++){
        while(!loGood){
            if(printSOUT)
                std::cout << "Simulating loci # " <<  i + 1 << std::endl;
            loGood = bdsaBDSim();
        }
        for(int j = 0; j < numGenes; j++){
            while(!gGood){
                if(printSOUT)
                    std::cout << "Simulating gene # " <<  j + 1 << " of loci # " << i + 1 << std::endl;
                gGood = coalescentSim();
            }
            geneTrees[i].push_back(geneTree);

            gGood = false;
        }
        locusTrees.push_back(lociTree);
        loGood = false;
    }
    return gGood;
}


/**
 * Calculates the tree depth as total branch length from root to tip
 * @return The tree depth as a double
 */
double Simulator::calcSpeciesTreeDepth(){
    return spTree->getTreeDepth();
}


/**
 * Calculates the tree depth of the tree with extinct tips pruned off
 * @return The tree depth of the pruned tree
 */
double Simulator::calcExtantSpeciesTreeDepth(){
    auto *tt = new SpeciesTree(rando, numTaxaToSim);
    spTree->getRootFromFlags(false);
    tt->setRoot(spTree->getExtantRoot());
    tt->setExtantRoot(tt->getRoot());
    tt->reconstructTreeFromSim(tt->getExtantRoot());
    double extTreeDepth = tt->getTreeDepth();
    delete tt;
    return extTreeDepth;
}

/**
 * Calculates the tree depth of the locus tree at index i in the vector of constructor LocusTree
 * @param i index of the locus tree
 * @return tree depth of the locus tree
 */
double Simulator::calcLocusTreeDepth(int i){
    return locusTrees[i]->getTreeDepth();
}

/**
 * Identifies the average number of transfers of all locus trees found in the vector of constructor LocusTree
 * @return The average transfers of the locus trees simulated for Simulator class
 */
int Simulator::findNumberTransfers(){
    auto numberTranfers = 0;
    for(auto & locusTree : locusTrees){
        numberTranfers += locusTree->getNumberTransfers();
    }
    return numberTranfers / (int) locusTrees.size();
} 
/**
 * Finds the average number of duplications in the simulated LocusTree for Simulator class
 * @return Average number of duplications
 */
int Simulator::findNumberDuplications(){
    auto numberDuplications = 0;
    for(auto & locusTree : locusTrees){
        numberDuplications += locusTree->getNumberDuplications();
    }
    return numberDuplications / (int) locusTrees.size();
}
/**
 * Finds the average number of losses in the simulated LocusTree for Simulator class
 * @return Average number of losses
 */

int Simulator::findNumberLosses(){
    auto numberLosses = 0;
    for(auto & locusTree : locusTrees){
        numberLosses += locusTree->getNumberLosses();
    }
    return numberLosses / (int) locusTrees.size();
}
/**
 * Finds average number of generations in gene trees for each locus tree
 * @return Average number of generations
 */
std::vector<double> Simulator::findNumberGenerations(){
    std::vector<double> numberGenerations;
    for(int i = 0; i < locusTrees.size(); i++){
        numberGenerations.push_back(0.0);
        for(int j = 0; j < geneTrees.size(); j++){
            if(!(geneTrees[i].empty())){
                numberGenerations[i] += geneTrees[i][j]->getTreeDepth() * popSize;
            }
        }
        numberGenerations[i] /= geneTrees.size();
    }
    return numberGenerations;
}
