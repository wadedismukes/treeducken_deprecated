//
//  main.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/7/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include <iostream>
#include "SpeciesTree.h"
#include "Simulator.h"
#include "Engine.h"

void printHelp();
void printSettings();

void printHelp(){
    std::cout << "\tHere are the available options that you can change (default values are in []):\n";
    std::cout << "\t\t-i    : input settings file \n";
    std::cout << "\t\t-o    : output file name prefix [= ""]\n";
    std::cout << "\t\t-nt   : number of extant taxa [= 100]\n";
    std::cout << "\t\t-r    : number of replicates [= 10]\n";
    std::cout << "\t\t-nloc : number of loci to simulate [= 0]\n";
    std::cout << "\t\t-sc   : tree scale [= 1.0]\n";
    std::cout << "\t\t-sd1  : seed 1 (use this if you only pass in one seed) \n";
    std::cout << "\t\t-sd1  : seed 2 \n";
    std::cout << "\t\t-sbr  : species birth rate [= 0.5]\n";
    std::cout << "\t\t-sdr  : species death rate [= 0.2]\n";
    std::cout << "\t\t-gbr  : gene birth rate [= 0.0]\n";
    std::cout << "\t\t-gdr  : gene death rate [= 0.0]\n";
    std::cout << "\t\t-lgtr : gene transfer rate [= 0.0]\n";
    std::cout << "\t\t-ipp  : individuals to sample per locus [= 0]\n";
    std::cout << "\t\t-ne   : effective population size per locus [= 0] \n";
    
}

void printSettings(std::string of, int nt, int r, int nloc, int ts, double sbr, double sdr,
                   double gbr, double gdr, double lgtr, int ipp, int ne){
    std::cout << "\t\toutput file name prefix         = " << of << "\n";
    std::cout << "\t\tNumber of extant taxa           = " << nt << "\n";
    std::cout << "\t\tNumber of replicates            = " << r << "\n";
    std::cout << "\t\tNumber of loci to simulate      = " << nloc << "\n";
    std::cout << "\t\tSpecies birth rate              = " << sbr << "\n";
    std::cout << "\t\tSpecies death rate              = " << sdr << "\n";
    std::cout << "\t\tGene birth rate                 = " << gbr << "\n";
    std::cout << "\t\tGene death rate                 = " << gdr << "\n";
    std::cout << "\t\tGene transfer rate              = " << lgtr << "\n";
    std::cout << "\t\tIndividuals to sample per locus = " << ipp << "\n";
    std::cout << "\t\tEffective pop size per locus    = " << ne << "\n";
}


int main(int argc, char * argv[]) {
    std::string setFileName;
    int mt;
    Engine *phyEngine;
    if(argc  == 1){
        printHelp();
        return 0;
    }
    else{
        std::string outName = "";
        int nt = 100, r = 10, nloc = 10, ipp = 0, ne = 0, sd1 = 0, sd2 = 0;
        double sbr = 0.5, sdr = 0.2, gbr = 0.0, gdr = 0.0, lgtr = 0.0, ts = 1.0;
        for (int i = 0; i < argc; i++){
                char *curArg = argv[i];
                if(strlen(curArg) > 1 && curArg[0] == '-'){
                    if(!strcmp(curArg, "-i")){
                        setFileName = argv[i+1];
                        std::ifstream settings(setFileName);
                        std::string line;
                        std::string comment = "#";
                        if(settings.is_open()){
                            while( getline (settings, line) ){
                                if(line.substr(0,1) != comment){
                                    if(line.substr(0,4) == "-sbr")
                                        sbr = atof(line.substr(5, std::string::npos - 1).c_str());
                                    else if(line.substr(0,4) == "-sdr")
                                        sdr = atof(line.substr(5, std::string::npos - 1).c_str());
                                    else if(line.substr(0,4) == "-gbr")
                                        gbr = atof(line.substr(5, std::string::npos - 1).c_str());
                                    else if(line.substr(0,4) == "-gdr")
                                        gdr = atof(line.substr(5, std::string::npos - 1).c_str());
                                    else if(line.substr(0,5) == "-lgtr")
                                        lgtr = atof(line.substr(6, std::string::npos - 1).c_str());
                                    else if(line.substr(0,3) == "-ne")
                                        ne = atof(line.substr(4, std::string::npos - 1).c_str());
                                    else if(line.substr(0,4) == "-ipp")
                                        ipp = atof(line.substr(5, std::string::npos - 1).c_str());
                                    else if(line.substr(0,3) == "-nt")
                                        nt = atoi(line.substr(4, std::string::npos - 1).c_str());
                                    else if(line.substr(0,3) == "-sc")
                                        ts = atof(line.substr(4, std::string::npos - 1).c_str());
                                    else if(line.substr(0,3) == "-nl")
                                        nloc = atof(line.substr(4, std::string::npos - 1).c_str());
                                    else if(line.substr(0,2) == "-r")
                                        r = atoi(line.substr(3, std::string::npos - 1).c_str());
                                    else if(line.substr(0,2) == "-o")
                                        outName = line.substr(3, std::string::npos - 1);
                                    else if(line.substr(0,3) == "-s1")
                                        sd1 = atoi(line.substr(4, std::string::npos - 1).c_str());
                                    else if(line.substr(0,3) == "-s2")
                                        sd2 = atoi(line.substr(4, std::string::npos - 1).c_str());
                                }
                            }
                            
                        }
                        else{
                            std::cerr << "The input file was unable to be opened.\n" << std::endl;
                            return 0;
                        }
                    }
                    else if(!strcmp(curArg, "-sbr"))
                        sbr = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-sdr"))
                        sdr = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-gbr"))
                        gbr = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-gdr"))
                        gdr = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-lgtr"))
                        lgtr = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-ne"))
                        ne = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-ipp"))
                        ipp = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-nt"))
                        nt = atoi(argv[i+1]);
                    else if(!strcmp(curArg, "-sc"))
                        ts = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-nl"))
                        nloc = atof(argv[i+1]);
                    else if(!strcmp(curArg, "-r"))
                        r = atoi(argv[i+1]);
                    else if(!strcmp(curArg, "-o"))
                        outName = argv[i+1];
                    else if(!strcmp(curArg, "-h")){
                        printHelp();
                        return 0;
                    }
                    else if(!strcmp(curArg, "-sd1"))
                        sd1 = atoi(argv[i+1]);
                    else if(!strcmp(curArg, "-sd2"))
                        sd2 = atoi(argv[i+1]);
                    else {
                        std::cout << "\n############################ !!! ###########################\n";
                        std::cout << "\n\n\tPerhaps you mis-typed something, here are the \n\tavailable options:\n";
                        printHelp();
                        std::cout << "\n############################ !!! ###########################\n";
                        return 0;
                    }
                }
        }

        if(sbr <= 0.0){
            std::cerr << "species birth rate of 0.0 is invalid, exiting...\n";
            return 0;
        }
        if(nloc > 0){
            if (gbr <= 0.0){
                if(gbr < 0.0){
                    mt = 1;
                    std::cout << "gene birth rate is a negative number, no loci or gene trees will be simulated.\n";
                }
                else{
                    if(ne >  0 && ipp > 0 && ipp <= ne){
                        mt = 3;
                        std::cout << "gene birth rate is 0.0, locus trees will match species trees.\n";
                    }
                    else{
                        std::cout << "gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                        std::cout << "population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                        mt = 2;
                    }
                }
                printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne);

            }
            else if (ne <= 0 || ipp <= 0 || ipp > ne){
                mt = 2;
                std::cout << "gene tree parameters are incorrectly specified.\n";
                std::cout << "population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne);

            }
            else{
                mt = 3;
                std::cout << "Simulating sets of three trees.\n";
            }
        }
        else{
            std::cout << "Number of loci to simulate is set to 0." << std::endl;
            mt = 1;
        }
        phyEngine = new Engine(outName,
                               mt,
                               sbr,
                               sdr,
                               gbr,
                               gdr,
                               lgtr,
                               ipp,
                               ne,
                               1.0,
                               sd1,
                               sd2,
                               ts,
                               r,
                               nt,
                               nloc);
        phyEngine->doRunRun();


    }
    
    delete phyEngine;
//    MbRandom rand;
//    int sd1 = 0;
//    int sd2 = 0;
//    if(sd1 > 0 && sd2 > 0)
//        rand.setSeed(sd1, sd2);
//    else
//        rand.setSeed();
//    seedType gs1, gs2;
//    rand.getSeed(gs1, gs2);
//    std::cout << "\nSeeds {" << gs1 << ", " << gs2 << "}" << std::endl;
//    
//    unsigned numTaxa = 4;
//    std::string sptree;
//    std::string locusTree;
//    std::string geneTree;
//    double br = 0.01;
//    double dr = 0.005;
//    double gbr = 0.02;
//    double gdr = 0.01;
//    double lgtr = 0.02;
//    unsigned ipp = 4;
//    unsigned popSize = 10000;
//    double genTime = 1;
//    unsigned numLoci = 10;
//    
//    Simulator *sim = new Simulator(&rand, numTaxa, br, dr, 1.0, numLoci, gbr, gdr, lgtr, ipp, popSize, genTime);
//    sim->simThreeTree();
//    sptree = sim->printSpeciesTreeNewick();
//    std::cout << "species tree: " << std::endl << sptree << std::endl <<  std::endl;
//    for(int i = 0; i < numLoci; i++){
//        locusTree = sim->printLocusTreeNewick(i);
//        std::cout << "Locus Tree: " << std::endl << locusTree << std::endl << std::endl;
//        geneTree = sim->printGeneTreeNewick(i);
//        std::cout << "Gene Tree: " << std::endl << geneTree << std::endl << std::endl;
//    }
    return 0;
}