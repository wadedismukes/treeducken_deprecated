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
    inputSpTree = "";
    simType = mt;
    spBirthRate = sbr;
    spDeathRate = sdr;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    doScaleTree = false;
    // treescale = 1.0;
    seedset = 0;

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
    for(std::vector<TreeInfo*>::iterator p=simSpeciesTrees.begin(); p != simSpeciesTrees.end(); ++p){
        delete (*p);
    }
    simSpeciesTrees.clear();
}


void Engine::doRunRun(){
    // double TS = 0.0;
    TreeInfo *ti = nullptr;
    for(int k = 0; k < numSpeciesTrees; k++){
        
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
  
        std::cout << "Simulating species tree replicate # " << k + 1 << std::endl;
        
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
        
        ti =  new TreeInfo(k, numLoci);
        ti->setWholeTreeStringInfo(treesim->printSpeciesTreeNewick());
        for(int i = 0; i < numLoci; i++){
            ti->setLocusTreeByIndx(k, treesim->printLocusTreeNewick(i));
            if(simType == 3){
                for(int j = 0; j < numGenes; j++){
                    ti->setGeneTreeByIndx(i, j, treesim->printGeneTreeNewick(i, j));
                    ti->setExtantGeneTreeByIndx(i, j, treesim->printExtantGeneTreeNewick(i, j));
                }
            }
        }
        simSpeciesTrees.push_back(ti);
        delete treesim;
        treesim = nullptr;
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
                // for(int j = 0; j < numGenes; j++){
                //     (*p)->writeGeneTreeFileInfoByIndx(d, i, j, outfilename);
                // }
                (*p)->writeGeneTreeFileInfo(d, i, numGenes, outfilename);
        }
    }
}


TreeInfo* Engine::findTreeByIndx(int i){
    TreeInfo *tf = 0;
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

SpeciesTree* Engine::buildTreeFromNewick(std::string spTreeStr){
    
}

void Engine::doRunSpTreeSet(){

    std::cout << "Setting species tree to this newick tree: " << inputSpTree << std::endl;

    TreeInfo *ti = nullptr;    
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
    

    treesim->setSpeciesTree(this->buildTreeFromNewick(inputSpTree));
    treesim->simLocusGeneTrees();


    ti =  new TreeInfo(0, numLoci);
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
    delete treesim;
    treesim = nullptr;

    this->writeTreeFiles();

}


/*
    TreeInfo functions to write tree information to file in various file formats
                                                                */
TreeInfo::TreeInfo(int idx, int nl){
    geneTrees.resize(nl);
    extGeneTrees.resize(nl);
    spTreeLength = 0.0;
    spTreeDepth = 0.0;
    spTreeNess = 0.0;
    spAveTipLen = 0.0;
    loTreeLength = 0.0;
    loTreeDepth = 0.0;
    loTreeNess = 0.0;
    loAveTipLen = 0.0;
    aveTMRCAGeneTree = 0.0;
}

TreeInfo::~TreeInfo(){
    gsaTrees.clear();
    locusTrees.clear();
    geneTrees.clear();
    extGeneTrees.clear();
    speciesTree.clear();
}

void TreeInfo::writeWholeTreeFileInfo(int spIndx, std::string ofp){
    std::string path = "";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    //path += tn.str();
    // path +=  "/";
    
    
    fn += "_" + tn.str() + ".sp.tre";
    path += fn;
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree wholeT_" << spIndx << " = ";
    out << getWholeSpeciesTree() << "\n";
    out << "end;";

}

void TreeInfo::writeLocusTreeFileInfoByIndx(int spIndx, int indx, std::string ofp){
    std::string path = "";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    // path += tn.str();
    // path += "/";
    
    
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
    std::string path = "";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    //path += tn.str();
   // path += "/";


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

void TreeInfo::writeExtGeneTreeFileInfo(int spIndx, int Lindx, int numGenes, std::string ofp){
    std::string path = "";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    //path += tn.str();
    // path += "/";
    fn += "genetrees_";

    fn += tn.str();
    
    tn.clear();
    tn.str(std::string());
 
    tn << Lindx;
    
    fn += "_" + tn.str() + ".tre";;
    tn.clear();
    tn.str(std::string());
    
    // tn << i;
    // fn += "_" + tn.str() + ".gen.tre";
    path += fn;
    
    std::ofstream out(path);
    // out << "#NEXUS\nbegin trees;\n";
    // for(int i = 0; i < numGenes; i++){
    //         out << "tree geneT_" << indx << " = ";
    //         out << getGeneTreeByIndx(Lindx, indx) << "\n";
    // }
    // out << "end;";
    
    
    out << "#NEXUS\nbegin trees;\n";
    for(int i = 0; i < numGenes; i++){
        out << "tree extGeneT_" << i << " = ";
        out << getExtGeneTreeByIndx(Lindx, i) << "\n";
    }
    out << "end;";

}

void TreeInfo::writeGeneTreeFileInfo(int spIndx, int Lindx, int numGenes, std::string ofp){
    std::string path = "";
    
    std::string fn = ofp;
    std::stringstream tn;
    
    tn << spIndx;
    //path += tn.str();
    // path += "/";
    fn += "genetrees_";

    fn += tn.str();
    
    tn.clear();
    tn.str(std::string());
 
    tn << Lindx;
    
    fn += "_" + tn.str() + ".tre";;
    tn.clear();
    tn.str(std::string());
    
    // tn << i;
    // fn += "_" + tn.str() + ".gen.tre";
    path += fn;
    
    std::ofstream out(path);
    // out << "#NEXUS\nbegin trees;\n";
    // for(int i = 0; i < numGenes; i++){
    //         out << "tree geneT_" << indx << " = ";
    //         out << getGeneTreeByIndx(Lindx, indx) << "\n";
    // }
    // out << "end;";
    
    
    out << "#NEXUS\nbegin trees;\n";
    for(int i = 0; i < numGenes; i++){
        out << "tree geneT_" << i << " = ";
        out << getGeneTreeByIndx(Lindx, i) << "\n";
    }
    out << "end;";

}
