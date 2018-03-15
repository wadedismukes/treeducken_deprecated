#!/bin/bash
# This file will read in settings from a simSettings.txt file, create a directory `sim_files` to store 
# a directory for each species tree to be simulated (output in Newick format)
# within each species tree directory a locus tree directory is created for each
# simulated locus tree and its corresponding gene tree to be stored in Newick 
# format. 

# Directions for writing the settings file can be found in README.md

echo "multiTree v. 0.1\n"

makedirs(){
  [ "$1" -gt 0 ] || return
  mkdir -v "speciestree_$1"
  makedirs $(( $1 - 1 ))
}

if [ "$#" -eq 3 ]; then
	[ -d sim_files ] || mkdir sim_files
	cd sim_files
	makedirs $2
  cd ..
	multiTree -i $1 -r $2
else
	echo "Incorrect number of parameters in config_sim.sh"
	echo "Correct usage is 'sh config_sim.sh simSettings.txt <#SpeciesTrees> <#LocusTrees>"
fi

