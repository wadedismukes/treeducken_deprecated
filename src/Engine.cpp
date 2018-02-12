//
//  Engine.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 1/28/18.
//  Copyright © 2018 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Engine.h"

Engine::Engine(std::string of, int mt, double sbr, double sdr, double gbr, double gdr, double lgtr, int ipp, int popsize, double genTime, int sd1, int sd2, double treescale, int reps, int ntax, int nloci){
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
    numTaxa = ntax;
    numLoci = nloci;
    if(sd1 > 0 && sd2 > 0)
        rando.setSeed(sd1, sd2);
    else
        rando.setSeed();
    seedType gs1, gs2;
    rando.getSeed(gs1, gs2);
    std::cout << "\nSeeds = {" << gs1 << ", " << gs2 << "}" << std::endl;
    
    if(treescale > 0.0)
        doScaleTree = true;
    else
        doScaleTree = false;
    
    if(spBirthRate <= 0.0){
        std::cerr << "ERROR: The starting birth rate is set to " << spBirthRate << " --  This is not right!" << std::endl;
        std::cerr << "No trees were generated." << std::endl;
        exit(1);

    }
    
}


Engine::~Engine(){

}


void Engine::doRunRun(){
    // to be written after designing an outfile scheme
}


void Engine::writeTreeFiles(){
    for(std::vector<TreeInfo *>::iterator p = simSpeciesTrees.begin(); p != simSpeciesTrees.end(); p++){
        (*p)->writeWholeTreeFileInfo(outfilename);
        for(int i = 0; i < numLoci; i++){
            (*p)->writeLocusTreeFileInfoByIndx(i, outfilename);
        }
        for(int j = 0; j < numGenes; j++){
            (*p)->writeGeneTreeFileInfoByIndx(j, outfilename);
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
void TreeInfo::writeWholeTreeFileInfo(std::string ofp){
    
}

void TreeInfo::writeLocusTreeFileInfoByIndx(int indx, std::string ofp){
    
}

void TreeInfo::writeGeneTreeFileInfoByIndx(int indx, std::string ofp){
    
}
