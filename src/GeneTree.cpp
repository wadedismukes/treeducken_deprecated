#include "GeneTree.h"
#include <iostream>

/**
 * @brief Constructor of the GeneTree class inherited from the Tree class
 *
 * @param MbRandom *p pointer to MrBayes RNG inherited from Tree class
 * @param nt Number of taxa as a positive integer inherited from Tree class
 * @param ipp Individuals per population to be sampled within each lineage branch
 * @param ne Population size as an integer
 * @param genTime Scaling factor for gene trees in units of generations per unit time
 */

GeneTree::GeneTree(MbRandom *p, unsigned nt, unsigned ipp, unsigned ne, double genTime) : Tree(p, nt){
    rando = p;
    numTaxa = nt;
    individualsPerPop = ipp;
    popSize = ne;
    generationTime = genTime;
    delete root;
}

/**
 * @brief Destructor of GeneTree class
 *
 *
 */

GeneTree::~GeneTree(){
    // for(std::vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); ++p){
    //     if((*p) != nullptr){
    //         delete (*p);
    //         (*p) = nullptr;
    //     }
    // }
    // clearNodes(extantRoot);
}


/**
 * @brief Function that takes vector of vector of indices of loci at the present and creates a vector of class Node at the start of GeneTree simulation
 *
 * @param extantLociInd Vector of a vector of indices to get the indices corresponding to extant loci
 * @param presentTime The time at present (or at end of locus tree)
 */

void GeneTree::initializeTree(std::vector< std::vector<int> > extantLociInd, double presentTime){
    for(auto & node : nodes){
        delete node;
        node = nullptr;
    }
    nodes.clear();
    for(auto & extantNode : extantNodes){
        delete extantNode;
        extantNode = nullptr;
    }
    extantNodes.clear();
    Node *p;
    int numberLociInPresent;
    numberLociInPresent = (int) extantLociInd[0].size();
    for(int i = 0; i < numberLociInPresent; i++){
        for(int j = 0; j < individualsPerPop; j++){
            p = new Node();
            p->setDeathTime(presentTime);
            p->setLindx(extantLociInd[0][i]);
            p->setIndx(extantLociInd[0][i]);
            p->setLdes(nullptr);
            p->setRdes(nullptr);
            p->setAnc(nullptr);
            p->setIsExtant(true);
            p->setIsTip(true);
            p->setIsExtinct(false);
            if(extantLociInd[0][i] == -1){
                this->setOutgroup(p);
                p->setName("OUT");
            }
            else{
                extantNodes.push_back(p);
            }
            nodes.push_back(p);
        }
    }
}

/**
 * @brief Gets the time to next coalescent event
 *
 *
 * @param n Number of currently living individuals in lineage
 * @return Time to event as a double
 */

double GeneTree::getCoalTime(int n){
    double ct;
    double lambda = (double)(n * (n - 1)) / (2 * popSize) ;
    ct = -log(rando->uniformRv()) / (lambda);
    return ct;
}

/**
 * @brief Function that coordinates the censored coalescent process
 * @details While time is not 0, gets times to events, checks which Nodes are alive randomly selects two and deletes them adding in a new Node
 *
 * @param startTime Left-bound on times based on censoring
 * @param stopTime Right-bound on times based on censoring
 * @param contempSpeciesIndx SpeciesIndx of lineage being coalesced in
 * @param ancSpIndx Index of the Ancestor of the Node of index contempSpecies Indx
 * @param chck True if last coalescent time, otherwise false
 * @return Bool indicating if number of nodes in Coalescent Process is now 1
 */

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
    for(auto it = extantNodes.begin(); it != extantNodes.end(); ++it){
        if((*it)->getLindx() == contempSpeciesIndx){
            extIndx = (int) std::distance(extantNodes.begin(), it);
            indInExtNodes.push_back(extIndx);
        }
    }
  //  std::cout << contempSpeciesIndx << "   ($)$)%    " << indInExtNodes.size() << std::endl;
    if(indInExtNodes.size() > 1){
        while(t > stopTime){
            t -= getCoalTime((int) indInExtNodes.size());
            // std::cout << t << std::endl;
            if(t < stopTime){
                {
                    allCoalesced = chck;
                    break;
                }
            }

            rightInd = rando->uniformRv(0, indInExtNodes.size() - 1);
            iter_swap(indInExtNodes.begin() + rightInd, indInExtNodes.begin());
            rightIndExtN = indInExtNodes[0];
            r = extantNodes[rightIndExtN];

            leftInd = rando->uniformRv(1, indInExtNodes.size() - 1);
            iter_swap(indInExtNodes.begin() + leftInd, indInExtNodes.begin() + 1);
            leftIndExtN = indInExtNodes[1];
            l = extantNodes[leftIndExtN];

            n = coalescentEvent(t, l, r);
            //iter_swap(extantNodes.begin() + rightIndExtN, extantNodes.end() - 1);
            //iter_swap(extantNodes.begin() + leftIndExtN, extantNodes.end() - 2);
            if(leftIndExtN > rightIndExtN){
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
            }
            else{
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
            }
            indInExtNodes.clear();
            //indInExtNodes.erase(indInExtNodes.begin(), indInExtNodes.begin() + 2);
            extantNodes.insert(extantNodes.begin(), n);

            // if(!(indInExtNodes.empty())){
            //     for(std::vector<int>::iterator it = indInExtNodes.begin(); it != indInExtNodes.end(); ++it){
            //         (*it) += 1;
            //     }
            // }
            for(auto it = extantNodes.begin(); it != extantNodes.end(); ++it){
                if((*it)->getLindx() == contempSpeciesIndx){
                    extIndx = (int) std::distance(extantNodes.begin(), it);
                    indInExtNodes.push_back(extIndx);
                }
            }
//            extantNodes.erase(extantNodes.end() - 2, extantNodes.end());
         //  std::cout << "size of extantNodes " << extantNodes.size() << std::endl;
            if(indInExtNodes.size() == 1){
                allCoalesced = true;
                break;
            }
        }
    }
    else if (indInExtNodes.size() == 1){
        allCoalesced = true;
        extantNodes[indInExtNodes[0]]->setLindx(ancSpIndx);
    }
    else{
        // std::cout << "this shouldn't even be happnening" << std::endl;
        allCoalesced = true;
    }

    if(allCoalesced){
        for(int indInExtNode : indInExtNodes){
          //  std::cout << "**************" << std::endl;
            extantNodes[indInExtNode]->setLindx(ancSpIndx);
        }
    }

    indInExtNodes.clear();
    return allCoalesced;
}

/**
 * @brief Function for coalescing of 2 nodes at time t to 1 node
 *
 *
 * @param t time of coalescent event
 * @param p pointer to left Node of coalescent event
 * @param q pointer to right Node of coalescent
 * @return pointer to new Node
 */

Node* GeneTree::coalescentEvent(double t, Node *p, Node *q){
    Node *n = new Node();
    n->setDeathTime(t);
    n->setLdes(p);
    n->setRdes(q);
    n->setIsExtant(false);
    n->setIsTip(false);
    n->setIsExtinct(false);
    n->setLindx(p->getLindx());
    n->setIndx(p->getIndex());
    nodes.push_back(n);

    p->setBirthTime(t);
    p->setAnc(n);
   // p->setSib(q);

    q->setBirthTime(t);
    q->setAnc(n);
    // q->setSib(p);


    return n;
}

/**
 * @brief rescales times based on timeMap
 *
 * @todo maybe delete this?
 * @param timeMap a multimap of indices with respective times
 * @return a multimap of recaled times with same indices
 */

std::multimap<int, double> GeneTree::rescaleTimes(const std::multimap<int, double>& timeMap){
    std::multimap<int, double> rescaledTimeMap;
    std::pair<int, double> p;
    for(auto & it : timeMap){
        p.first = it.first;
        p.second = (it.second);
        rescaledTimeMap.insert(p);
    }

    return rescaledTimeMap;

}

/**
 * @brief Coalescent process that occurs at the root
 * @details Loops through any remaining nodes and coalesces them at times given by getCoalTime. Note that this means negative branch lengths are possible - these indicate events occuring within the stem.
 *
 * @todo get rid of ogf because it is clunky and doesn't work?
 * @param startTime time of root (final left-bound)
 * @param ogf fraction of outgroup
 */

void GeneTree::rootCoalescentProcess(double startTime, double ogf){
    int leftInd, rightInd;
    Node *l, *r;
    Node *n;
    double t = startTime;
    // search extantNodes for members with Lindx = contempSpecisIndx
    std::vector<int> indInExtNodes;
    for(auto & extantNode : extantNodes){
        extantNode->setLindx(0);
    }
    while(extantNodes.size() > 1){
        t -= getCoalTime((int) extantNodes.size());

        rightInd = rando->uniformRv(0, extantNodes.size() - 1);
        r = extantNodes[rightInd];
        extantNodes.erase(extantNodes.begin() + rightInd);

        leftInd = rando->uniformRv(0, extantNodes.size() - 1);
        l = extantNodes[leftInd];
        extantNodes.erase(extantNodes.begin() + leftInd);

        n = coalescentEvent(t, l, r);
        extantNodes.push_back(n);
    }
    if(ogf == 0.0){
        extantNodes[0]->setAsRoot(true);
        this->setRoot(extantNodes[0]);
    }
    else{
        Node *nRoot = new Node();
        t -= getCoalTime(2);
        extantNodes[0]->setBirthTime(t);
        nRoot->setBirthTime(t);
        nRoot->setLdes(extantNodes[0]);
        nRoot->setRdes(this->getOutgroup());
        nRoot->setDeathTime(extantNodes[0]->getBirthTime());
        nRoot->setAsRoot(true);

        this->setRoot(nRoot);
        this->getOutgroup()->setBirthTime(extantNodes[0]->getBirthTime());
        this->getOutgroup()->setAnc(nRoot);
        extantNodes[0]->setAnc(nRoot);
        nodes.push_back(nRoot);

    }


    this->setBranchLengths();
}

/**
 * @brief Recursive function for rescaling times to remove negative branch lengths
 * @todo check where this is called
 *
 *
 * @param r Pointer to root Node (position one of nodes vector[Node*])
 * @param add The number to add to put the coalescent tree root at 0
 */

void GeneTree::recursiveRescaleTimes(Node* r, double add){
    if(r != nullptr){
        if( r->getRdes() == nullptr){
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

/**
 * @brief Sets branch lengths of the GeneTree
 * @details This is the same function as SpeciesTree and LocusTree
 */
void GeneTree::setBranchLengths(){
    double brlen;
    for(auto & node : nodes){
        brlen = node->getDeathTime() - node->getBirthTime();
        node->setBranchLength(std::abs(brlen));
    }
}

/**
 * @brief Adds in extinct species so that coalescent events can occur on these trees with same number of individuals as ipp indicates
 *
 *
 * @param bt Time of death of the lineage (birth in coalescent process)
 * @param indx index of the extinct lineage in nodes[*Node] vector of SpeciesTree class
 */

void GeneTree::addExtinctSpecies(double bt, int indx){
    Node *p;
    for(int i = 0; i < individualsPerPop; i++){
        p = new Node();
        p->setDeathTime(bt);
        p->setIndx(indx);
        p->setLindx(indx);
        p->setLdes(nullptr);
        p->setRdes(nullptr);
        p->setAnc(nullptr);
        p->setIsExtant(false);
        p->setIsTip(true);
        p->setIsExtinct(true);
        extantNodes.push_back(p);
        nodes.push_back(p);

    }
    //delete p;
}


/**
 * @brief Resets indices of LocusTree class to indices of SpeciesTree class
 *
 *
 * @param spToLocusMap map of species indices to locus indices
 */
void GeneTree::setIndicesBySpecies(std::map<int, int> spToLocusMap){
    int indx;
    int spIndx;
    for(auto & node : nodes){
        if(node->getIsTip()){
            indx = node->getIndex();
            spIndx = spToLocusMap.find(indx)->second;
            node->setIndx(spIndx);
        }
        else{
            indx = node->getLindx();
            spIndx = spToLocusMap.find(indx)->second;
            node->setIndx(spIndx);
        }
    }
    this->setTreeTipNames();
}

/**
 * @brief Printing function for GeneTree class
 * @details recursively writes tree by calling recGetNewickTree
 * @return Newick string of the tree
 */

std::string GeneTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string geneTreeString = ss.str();
    return geneTreeString;
}

/**
 * @brief Prints tree found in GeneTree class with only extant tips
 *
 * @todo redundant code; delete!
 * @return Newick string of the tree
 */
std::string GeneTree::printExtantNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string geneTreeString = ss.str();
    return geneTreeString;
}


void GeneTree::recGetExtNewickTree(Node *p, std::stringstream &ss){
    if(p->getRdes() == nullptr){
        if(p->getIsExtant())
            ss << p->getName();
    }
    else{
    // if(p != NULL){
        int flag = p->getFlag();
        if(flag == 2){
            ss << "(";
            recGetExtNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ",";
            recGetExtNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ")";
        }
        else{
            // if(p->getRdes() == NULL){
            //     ss << p->getName();
            // }
            // else{
            if(p->getLdes()->getIsExtinct()){
                recGetExtNewickTree(p->getRdes(), ss);
            }
            if(p->getRdes()->getIsExtinct()){
                recGetExtNewickTree(p->getLdes(), ss);
            }
            //}
        }
   }
}

/**
 * @brief Recursively writes a Newick string using std::stringstream based on tree structure in nodes
 *
 *
 * @param p Node at present to be written to ss
 * @param ss stringstream of the string to be printed out
 */
void GeneTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != nullptr){
        if( p->getRdes() == nullptr)
            ss << p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ")";
        }
    }
}

/**
 * @brief Loops through the vector of class Nodes and writes names to tips
 *
 */
void GeneTree::setTreeTipNames(){
    int indNumber = 0;
    std::stringstream tn;
    std::string name;
    int indx;
    for(auto & node : nodes){
        if(node->getIsTip()){
            indx  = node->getIndex();
            tn << indx;
            name = tn.str();
            tn.clear();
            tn.str(std::string());
            indNumber++;
            tn << indNumber;
            name += "_" + tn.str();
            tn.clear();
            tn.str(std::string());
            if(node == this->getOutgroup())
                node->setName("OUT");
            else
                node->setName(name);
            if(indNumber == individualsPerPop)
                indNumber = 0;
        }

    }

}
