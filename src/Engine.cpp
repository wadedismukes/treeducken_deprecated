//
//  Engine.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 1/28/18.
//  Copyright © 2018 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Engine.h"

Engine::Engine(std::string of, int mt, double sbr, double sdr, double gbr, double gdr, double lgtr, int ipp, int popsize, double genTime, int sd1, int sd2, double treescale, int reps, int ntax, int nloci, int ngen, double og){
    outfilename = of;
    simType = mt;
    spBirthRate = sbr;
    spDeathRate = sdr;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    individidualsPerPop = ipp;
    populationSize = popsize;
    generationTime = genTime;
    numSpeciesTrees = reps;
    proportionToSample = 1.0;
    numTaxa = ntax;
    numLoci = nloci;
    numGenes = ngen;
    outgroupFrac = og;
    if(sd1 > 0 && sd2 > 0)
        rando.setSeed(sd1, sd2);
    else
        rando.setSeed();
    seedType gs1, gs2;
    rando.getSeed(gs1, gs2);
    std::cout << "\nSeeds = {" << gs1 << ", " << gs2 << "}" << std::endl;
    
    
}


Engine::~Engine(){
    if(!(simSpeciesTrees.empty()))
        simSpeciesTrees.clear();
}


void Engine::doRunRun(){
    // to be written after designing an outfile scheme
    double TS = 0.0;
    for(int i = 0; i < numSpeciesTrees; i++){
        
        Simulator *treesim = new Simulator(&rando,
                                           numTaxa,
                                           spBirthRate,
                                           spDeathRate,
                                           1.0,
                                           numLoci,
                                           geneBirthRate,
                                           geneDeathRate,
                                           transferRate,
                                           individidualsPerPop,
                                           populationSize,
                                           generationTime,
                                           numGenes,
                                           outgroupFrac);
        std::cout << "simulating species tree replicate #" << i + 1 << std::endl;
        switch(simType){
            case 1:
                treesim->simSpeciesTree();
                break;
            case 2:
                treesim->simSpeciesLociTrees();
                break;
            case 3:
                treesim->simThreeTree();
                break;
            default:
                treesim->simSpeciesTree();
                break;
        }
        
        TreeInfo *ti = new TreeInfo(i, numLoci);
        ti->setWholeTreeStringInfo(treesim->printSpeciesTreeNewick());
        for(int i = 0; i < numLoci; i++){
            ti->setLocusTreeByIndx(i, treesim->printLocusTreeNewick(i));
            if(simType == 3){
                for(int j = 0; j < numGenes; j++){
                    ti->setGeneTreeByIndx(i, j, treesim->printGeneTreeNewick(i, j));
                    ti->setExtantGeneTreeByIndx(i, j, treesim->printExtantGeneTreeNewick(i, j));
                }
            }
        }
        simSpeciesTrees.push_back(ti);
        // TODO: rewrite the set whole tree string info to insert the string for each tree in a run
        // TODO: write function defintions for everything in engine
    }

    this->writeTreeFiles();
}


void Engine::writeTreeFiles(){
    
    for(std::vector<TreeInfo *>::iterator p = simSpeciesTrees.begin(); p != simSpeciesTrees.end(); p++){
        int d = (int) std::distance(simSpeciesTrees.begin(), p);
        (*p)->writeWholeTreeFileInfo(d, outfilename);
        for(int i = 0; i < numLoci; i++){
            (*p)->writeLocusTreeFileInfoByIndx(d, i, outfilename);
            if(simType == 3)
                for(int j = 0; j < numGenes; j++){
                    (*p)->writeGeneTreeFileInfoByIndx(d, i, j, outfilename);
                }
        }
    }
}


TreeInfo* Engine::findTreeByIndx(int i){
    TreeInfo *tf;
    int count = 0;
    for(std::vector<TreeInfo*>::iterator it = simSpeciesTrees.begin(); it != simSpeciesTrees.end(); ++it){
        if(count == i){
            tf = (*it);
            break;
        }
        else
            count++;
    }
    return tf;
}


void Engine::calcAverageRootAgeSpeciesTrees(){
    std::ofstream out;
    out.open("Average_root_depths_spTree.out");
    double sumRH = 0.0;
    for(std::vector<TreeInfo*>::iterator p = simSpeciesTrees.begin(); p != simSpeciesTrees.end(); ++p){
        sumRH += (*p)->getSpeciesTreeDepth();
        out << (*p)->getSpeciesTreeDepth() << "\n";
    }
    out.close();
    
    std::cout << "\n #### Average Root Age of Species Trees is " << (double) sumRH / numSpeciesTrees << " #####" << std::endl;
    
    
}


/*
    TreeInfo functions to write tree information to file in various file formats
                                                                                    */
void TreeInfo::writeWholeTreeFileInfo(int spIndx, std::string ofp){
    std::string path = "./sim_files/speciestree_";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    path += tn.str() + "/";
    
    
    fn += "_" + tn.str() + ".sp.tre";
    path += fn;
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree wholeT_" << spIndx << " = ";
    out << getWholeSpeciesTree() << "\n";
    out << "end;";

}

void TreeInfo::writeLocusTreeFileInfoByIndx(int spIndx, int indx, std::string ofp){
    std::string path = "./sim_files/speciestree_";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    path += tn.str();
    path += "/";
    
    
    fn += "_" + tn.str();
    
    tn.clear();
    tn.str(std::string());
    
    tn << indx;
    
    fn += "_" + tn.str() + ".loc.tre";
    path += fn;
    
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree locT_" << indx << " = ";
    out << getLocusTreeByIndx(indx) << "\n";
    out << "end;";
}

void TreeInfo::writeGeneTreeFileInfoByIndx(int spIndx, int Lindx, int indx, std::string ofp){
    std::string path = "./sim_files/speciestree_";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    path += tn.str();
    path += "/";


    fn += "_" + tn.str();
    
    tn.clear();
    tn.str(std::string());
 
    tn << Lindx;
    
    fn += "_" + tn.str();
    tn.clear();
    tn.str(std::string());
    
    tn << indx;
    fn += "_" + tn.str() + ".gen.tre";
    path += fn;
    
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree geneT_" << indx << " = ";
    out << getGeneTreeByIndx(Lindx, indx) << "\n";
    out << "tree extGeneT_" << indx << " = ";
    out << getExtGeneTreeByIndx(Lindx, indx) << "\n";
    out << "end;";

}
