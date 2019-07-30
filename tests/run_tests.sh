#!/bin/bash

cd test-1/
treeducken -i 025-turnover-settings.txt
mkdir -p output/
mv *.tre output/
cd ../test-2/

treeducken -i 0-5-trnsfr-settings.txt
mkdir -p output/
mv *.tre output/
cd ../test-3/

treeducken -i popsize-100-settings.txt
mkdir -p output/
mv *.tre output/
cd ..

# run R script to check files