#!/bin/bash
# This file will read in settings from a simSettings.txt file, create a directory `sim_files` to store 
# a directory for each species tree to be simulated (output in Newick format)
# within each species tree directory a locus tree directory is created for each
# simulated locus tree and its corresponding gene tree to be stored in Newick 
# format. 

# Directions for writing the settings file can be found in README.md


makedirs(){
  [ "$1" -gt -1 ] || return
  mkdir -v "speciestree_$1"
  makedirs $(( $1 - 1 ))
}

[ -d sim_files ] || mkdir sim_files
cd sim_files
makedirs $(( $1 - 1 ))

