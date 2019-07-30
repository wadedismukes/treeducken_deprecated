#ifndef LocusTree_h
#define LocusTree_h

#include "SpeciesTree.h"
#include <algorithm>
#include <set>
/**
 * @brief LocusTree class which is a child of the LocusTree class. 
 * @details LocusTree is produced via functions within the Simulator class. Within the tree structure of the SpeciesTree class
 *
 */
class LocusTree : public Tree 
{
    private:
        double stopTime;
        double geneBirthRate, geneDeathRate, transferRate;
        double currentTime;
        unsigned numTaxa;
        unsigned numTransfers;
        unsigned numDuplications;
        unsigned numLosses;

    public:
        LocusTree(MbRandom *rando, unsigned nt, double stop, double gbr, double gdr, double lgtr);
        virtual         ~LocusTree();
        virtual double  getTimeToNextEvent();
        virtual void    lineageBirthEvent(unsigned indx);
        virtual void    lineageDeathEvent(unsigned indx);
        virtual void    setNewLineageInfo(int indx, Node *r, Node *s);
        void    lineageTransferEvent(int indx);
        void    ermEvent(double ct);
    
        int     speciationEvent(int indx, double time, std::pair<int,int> sibs);
        void    extinctionEvent(int indx, double time);

        std::string   printNewickTree();
        void    setTreeTipNames();
        void    recTipNamer(Node *p, unsigned &copyNumber);
        static void    recGetNewickTree(Node *r, std::stringstream &ss);
        void    setBranchLengths();
        void    setPresentTime(double currentT);
        void    setStopTime(double st) {stopTime = st;}
        double  getCurrentTime() { return currentTime; }
        void    setCurrentTime(double ct) {currentTime = ct; }
        int     getNumberTransfers();
        int     getNumberDuplications();
        int     getNumberLosses();
        std::map<int,double>     getBirthTimesFromNodes();
        std::set<int>            getExtLociIndx();
        std::set<int>            getCoalBounds();
        std::multimap<int,double>     getDeathTimesFromNodes();
        std::multimap<int,double>     getDeathTimesFromExtinctNodes();
        std::map<int,int>             getLocusToSpeciesMap();
        std::vector< std::vector<int> >     getExtantLoci(const std::set<double, std::greater<double> >& epochSet);
        std::vector< std::string >    printSubTrees();
        int     postOrderTraversalStep(int indx);
    
    
    

};
#endif /* LocusTree_h*/
