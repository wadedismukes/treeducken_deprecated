#include <iostream>
#include "LocusTree.h"

/**
 * @brief Constructor for the locus tree class 
 * @details Inherited from the tree class, structure is based on the nodes vector of class Node* of the SpeciesTree class.
 * 
 * @param p MrBayes based RNG seed
 * @param nt Number of tips in species tree as unsigned int
 * @param stop time to stop at determined by SpeciesTree as double
 * @param gbr Gene birth rate as double
 * @param gdr Gene death rate as double
 * @param lgtrate Lateral gene transfer rate 
 */
LocusTree::LocusTree(MbRandom *p, unsigned nt, double stop, double gbr, double gdr, double lgtrate) : Tree(p, nt, 0.0){
    rando = p;
    numTaxa = 1;
    stopTime = stop;
    currentTime = 0.0;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtrate;
    numTransfers = 0;
    numLosses = 0;
    numDuplications = 0;
    getRoot()->setLindx(0);
}

/**
 * @brief Destructor of LocusTree class
 */

LocusTree::~LocusTree(){

}

/**
 * @brief Function that sets information of Node during simulation
 * @details for both left and right lineages sets information for descendants
 * 
 * @param indx Index of node that is splitting in the vector of class Node*
 * @param r right Node*
 * @param l left Node*
 */

void LocusTree::setNewLineageInfo(int indx, Node *r, Node *l){
    extantNodes[indx]->setLdes(l);
    extantNodes[indx]->setRdes(r);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsDuplication(true);

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


    extantNodes.push_back(r);
    extantNodes.push_back(l);
    r->setLindx((int)nodes.size());

    nodes.push_back(r);
    l->setLindx((int) nodes.size());
    nodes.push_back(l);
    extantNodes.erase(extantNodes.begin() + indx);

    numExtant = (int)extantNodes.size();
}

/**
 * @brief Function for birth event of Node at indx in nodes vector
 * 
 * 
 * @param indx Index of node being split into right and left descendant
 */

void LocusTree::lineageBirthEvent(unsigned indx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    setNewLineageInfo(indx, right, sis);
    numDuplications += 1;
}

/**
 * @brief Function for lineage extinction
 * 
 * 
 * @param indx Index of lineage that is going extinct
 */

void LocusTree::lineageDeathEvent(unsigned indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numLosses += 1;
    numExtant = (int) extantNodes.size();
}

/**
 * @brief Function for transfer of one locus lineage from one species lineage to another
 *
 * 
 * @param indx Index of locus lineage that is to be transferred
 */

void LocusTree::lineageTransferEvent(int indx){
    unsigned spIndxD = extantNodes[indx]->getIndex();
    int allSameSpInExtNodes = 0;
//    if(numTaxa < 5){
        for(std::vector<Node*>::iterator p = extantNodes.begin(); p != extantNodes.end(); ++p){
            unsigned checkInd = (*p)->getIndex();
            if(checkInd != spIndxD){
                allSameSpInExtNodes++;
            }
        }
        if(allSameSpInExtNodes == 0)
            return;
 //   }
    //first a birth event
    Node *donor, *rec;
    donor = new Node();
    rec = new Node();
    numTransfers++;
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


    //extantNodes[indx] `
    extantNodes[indx]->setLdes(rec);
    extantNodes[indx]->setRdes(donor);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setFlag(1);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsDuplication(true);

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
    // rec->setAnc(extantNodes[recIndx.first]->getAnc());
    rec->setSib(NULL);

    speciesIndx.clear();
    extantNodes[recIndx.first]->setLdes(NULL);
    extantNodes[recIndx.first]->setRdes(NULL);
    extantNodes[recIndx.first]->setDeathTime(currentTime);
    extantNodes[recIndx.first]->setFlag(1);
    extantNodes[recIndx.first]->setIsExtant(false);
    extantNodes[recIndx.first]->setIsExtinct(true);
    extantNodes[recIndx.first]->setIsTip(true);
    iter_swap(extantNodes.begin() + recIndx.first, extantNodes.end()-1);
    iter_swap(extantNodes.begin() + indx, extantNodes.end()-2);

    extantNodes.erase(extantNodes.end()-2, extantNodes.end());
    extantNodes.push_back(rec);
    extantNodes.push_back(donor);

    rec->setLindx((int)nodes.size());
    nodes.push_back(rec);

    donor->setLindx((int)nodes.size());
    nodes.push_back(donor);

    numExtant = (int) extantNodes.size();
}

/**
 * @brief Function to determine time till the next event based on rates of events
 * 
 * @return Double of time till the next event
 */
double LocusTree::getTimeToNextEvent(){
    double sumrt = geneBirthRate + geneDeathRate + transferRate;
    double returnTime = 0.0;
    returnTime = -log(rando->uniformRv()) / (sumrt);
    currentTime += returnTime;
    return returnTime;
}

/**
 * @brief Event function to determine which event based on rates
 * 
 * 
 * @param ct Time as double of event
 */
void LocusTree::ermEvent(double ct){
    double relBr = geneBirthRate / (geneDeathRate + geneBirthRate + transferRate);
    double relLGTr = transferRate / (geneBirthRate + geneDeathRate + transferRate) + relBr;
    double whichEvent = rando->uniformRv();
    unsigned long extantSize = extantNodes.size();
    unsigned nodeInd = rando->uniformRv(0, extantSize - 1);
    currentTime = ct;
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


/**
 * @brief Function that splits ALL locus lineages in a species lineage into species descendants
 * @details Based on time splits all locus lineages in species lineage at indx and creates new lineages to be copied in descedants given in sibs
 * 
 * @param indx Index of species that speciates at time time
 * @param time Time speciation event is occurring
 * @param sibs Pair of indices of descendant lineages of species lineage at indx
 * @return Number of lineages out
 */
int LocusTree::speciationEvent(int indx, double time, std::pair<int,int> sibs){
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
            r->setIndx(sibs.second);

            l->setLdes(NULL);
            l->setRdes(NULL);
            l->setSib(r);
            l->setAnc((*it));
            l->setBirthTime(time);
            l->setIsTip(true);
            l->setIsExtinct(false);
            l->setIsExtant(true);
            l->setIndx(sibs.first);

            (*it)->setLdes(l);
            (*it)->setRdes(r);
            (*it)->setDeathTime(time);
            (*it)->setIsTip(false);
            (*it)->setIsExtant(false);
            r->setLindx((int)nodes.size());
            nodes.push_back(r);
            l->setLindx((int)nodes.size());
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

/**
 * @brief Extinction event for locus lineages in species lineage
 * @details For extinction events in SpeciesTree at time time for lineage at indx deletes all locus lineages
 * 
 * @param indx Index of species lineage that is going extinct
 * @param time double giving the time of extinction
 */
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
            ++it;
        }
    }
    numTaxa--;
}

/**
 * @brief Sets new indices of locus lineages after speciation events
 *
 * 
 * @param indx Index of species lineage that speciates
 * @param sibs  Pair of indices giving left as 0 and right as 1 
 * @param count Number of locus lineages to be set 
 */

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
            ++it;

    }
}

/**
 * @brief Recursive function for printing Newick tree
 * 
 * 
 * @param p Node to have index and branch length printed
 * @param ss Stringstream where the Newick string tree is written to 
 */

void LocusTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL)
            ss << p->getName();
        else{
            if(p->getFlag() == 1){
                ss << "(";
                ss << "(";
                recGetNewickTree(p->getLdes(), ss);
                ss << "[&index=" << p->getLdes()->getIndex() << "]" <<  ":" << p->getLdes()->getBranchLength();
                ss << ",";
                recGetNewickTree(p->getRdes(), ss);
                ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
                ss << ")";
                ss << p->getName() << ":" << "0.0";
                ss << ")";
            }
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
}

/**
 * @brief Sets final time at end of locus tree simulation
 * 
 * @param currentT Current time as a double in simulation 
 */

void LocusTree::setPresentTime(double currentT){
    for(std::vector<Node*>::iterator it = nodes.begin(); it !=  nodes.end(); ++it){
        if((*it)->getIsExtant())
            (*it)->setDeathTime(currentT);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}

//TODO: definitely make this a thing
std::vector<std::string> LocusTree::printSubTrees(){
    std::vector<std::string> subTrees;
    std::stringstream ss;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsDuplication()){
            recGetNewickTree((*it),ss);
            ss << ";";
            subTrees.push_back(ss.str());
            ss.clear();
            ss.str(std::string());
        }
    }
    return subTrees;
}

/**
 * @brief Function to call recGetNewickTree 
 * @return Returns Newick string of a tree 
 */
std::string LocusTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string loTree = ss.str();
    return loTree;
}

/**
 * @brief Sets the tip names of a tree according to species indices
 * 
 */

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
                    tn << (*it)->getRdes()->getIndex();
                    copyNumberCounts.push_back((*it)->getIndex());
                    name = "TR" + tn.str();
                    tn.clear();
                    tn.str(std::string());
                    tn << (*it)->getLdes()->getIndex();
                    name += "->" + tn.str();
                    (*it)->setName(name);
                    tn.clear();
                    tn.str(std::string());

                }

            }
        }
    }
}

// TODO: this isn't being used and might need to be deleted
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


/**
 * @brief Sets branch lengths for a LocusTree
 * @details Same for the other classes inherited from Tree
 */

void LocusTree::setBranchLengths(){
    double brlen;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        brlen = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(brlen);
    }
}

/**
 * @brief Function to get death times of all LocusTree nodes
 * @return Multimap containing pairs of species tree nodes with death times
 */
std::multimap<int, double> LocusTree::getDeathTimesFromNodes(){
    int locusIndx;
    double deathTime;
    std::multimap<int,double> deathTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        locusIndx = (int)(*it)->getLindx();
        deathTime = (*it)->getDeathTime();
        deathTimeMap.insert(std::pair<int,double>(locusIndx, deathTime));
    }
    return deathTimeMap;
}

/**
 * @brief Function to get death times of LocusTree nodes that are going extinct
 * @return Multimap containing pairs of species tree nodes with death times
 */

std::multimap<int, double> LocusTree::getDeathTimesFromExtinctNodes(){
    int locusIndx;
    double deathTime;
    std::multimap<int, double> deathTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsExtinct()){
            locusIndx = (int)(*it)->getLindx();
            deathTime = (*it)->getDeathTime();
            deathTimeMap.insert(std::pair<int,double>(locusIndx, deathTime));
        }

    }

    return deathTimeMap;
}

/**
 * @brief Function to get birth times of all LocusTree nodes
 * @return Multimap containing pairs of species tree nodes with birth times
 */
std::map<int, double> LocusTree::getBirthTimesFromNodes(){
    int locusIndx;
    double birthTime;
    std::map<int,double> birthTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        locusIndx = (int)(*it)->getLindx();
        birthTime = (*it)->getBirthTime();
        birthTimeMap.insert(std::pair<int,double>(locusIndx, birthTime));
    }
    if(this->getOutgroup() != NULL)
        birthTimeMap.insert(std::pair<int,double> (-1, this->getOutgroup()->getBirthTime()));
    return birthTimeMap;
}

/**
 * @brief Function to retrieve a vector of which loci are alive at which times slices
 * @details Returns a vector of vector. The outer vector is of length equal to the number of 
 epochs" defined by either branching events in locus tree or death events. The inner vectors are the indices alive within each epoch
 * 
 * @param epochs Set of epochs sorted from greatest to least
 * @return Vector of epochs with members vectors of indices
 */)
std::vector< std::vector<int> > LocusTree::getExtantLoci(std::set<double, std::greater<double> > epochs){

    int locusIndx;
    int epCount = 0;
    int numEpochs = (int) epochs.size();
    std::vector< std::vector<int> > locusInEpoch(numEpochs);
    for(std::set<double, std::greater<double> >::iterator epIt = epochs.begin(); epIt != epochs.end(); ++epIt){
        for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
            if(epCount == 0){
                if((*it)->getIsExtant()){
                    locusIndx = (*it)->getLindx();
                    locusInEpoch[epCount].push_back(locusIndx);

                }
            }
            else{
                if((*it)->getDeathTime() >= (*epIt)){
                    locusIndx = (*it)->getLindx();
                    locusInEpoch[epCount].push_back(locusIndx);
                }
            }
        }



        epCount++;
    }
    if(this->getOutgroup() != NULL){
        locusInEpoch[0].push_back(-1);
    }
    return locusInEpoch;
}

/**
 * @brief Function to go in postorder (rootward) traversal
 * 
 * 
 * @param indx Starting index
 * @return Index of ancestor of indx
 */
int LocusTree::postOrderTraversalStep(int indx){
    Node* anc;
    int ancIndx;
    anc = nodes[indx]->getAnc();
    if(anc != NULL)
        ancIndx = anc->getLindx();
    else
        ancIndx = 0;


    return ancIndx;
}

/**
 * @brief Function to get a map of loci indices associated with species indices
 * 
 * @return A map of loci indices with species indices
 */

std::map<int,int> LocusTree::getLocusToSpeciesMap(){
    std::map<int,int> locusToSpecies;
    int spID, loID;
    std::pair<int,int> pp;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        loID = (*it)->getLindx();
        spID = (*it)->getIndex();
        pp.first = loID;
        pp.second = spID;
        // std::cout << "locus id " << pp.first << " species id " << pp.second << std::endl;
        locusToSpecies.insert(locusToSpecies.end(), pp);
    }
    return locusToSpecies;
}

/**
 * @brief Gets indices of lineages that go extinct before present
 * @return Set of indices of Node* classes from LocusTree that are Extinct and Tips
 */
std::set<int> LocusTree::getExtLociIndx(){
    std::set<int> doomedLoci;
    int indx;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsExtinct() && (*it)->getIsTip()){
            indx = (int) std::distance(nodes.begin(), it);
            doomedLoci.insert(doomedLoci.begin(), indx);
        }
    }
    return doomedLoci;
}


std::set<int> LocusTree::getCoalBounds(){
    std::set<int> coalBoundLoci;
    int indx;
    for(std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsDuplication()){
            indx = (int) std::distance(nodes.begin(), it);
            coalBoundLoci.insert(coalBoundLoci.begin(), indx);
        }
    }
    return coalBoundLoci;
}
/**
 * @brief Returns the number of transfers that have occured in LocusTree history
 * @return Integer of number of transfers
 */
int LocusTree::getNumberTransfers(){
    return numTransfers;
}
/**
 * @brief Returns the number of gene duplications (births) that have occured in LocusTree history
 * @return Integer of number of duplications
 */
int LocusTree::getNumberDuplications(){
    return numDuplications;
}
/**
 * @brief Returns the number of gene losses (deaths) that have occured in LocusTree history
 * @return Integer of number of gene losses
 */
int LocusTree::getNumberLosses(){
    return numLosses;
}
