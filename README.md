# `treeducken` 
`treeducken`, named after the infamous [dish](https://en.wikipedia.org/wiki/Turducken), is for the simulation of phylogenies under the three tree model (Rasmussen &amp; Kellis 2012). These three trees represent the species tree, locus (i.e. gene family) tree, and the gene tree. If you prefer, you may think of the species tree as the turkey, locus tree as the duck, and gene tree as the chicken in our phylogenetic three-bird roast.

--------

# Brief description of the simulator

The simulations require the input of the following parameters:

* number of species to simulate
* species birth rate
* species death rate
* gene birth rate
* gene death rate
* gene transfer rate
* individuals to sample per population
* effective population size
* number of loci to simulate on the species tree
* number of genes to simulate for each locus

After inputting these parameters the simulation proceeds as follows. First, a species tree is simulated under the constant rate birth-death process until the specified number of species is obtained using the generalized sampling algorithm given in Hartmann 2010. Next, a set of locus tree is forward simulated along the species tree using a birth-death process where transfers are treated as birth events where one loci is transferred to a randomly-selected contemporaneous species, and the other loci stays within the donor lineage. Note that this does simulate along lineages that are bound for extinction. The locus tree simulation process is repeated for the number of loci to simulate on the species tree. For each of those simulated locus trees a gene tree is backward simulated using the multi-locus coalescent. For this multilocus coalescent each tip in the tree starts with the individuals to sample per population. 


# Installing `treeducken`

To install `treeducken`, clone the repository by doing the following:

```
	git clone https://github.com/wadedismukes/treeducken.git
```

Then change directory to `treeducken/src` and run the following command:

```
	make install
```

`treeducken` is now installed in the `treeducken` directory. To uninstall simply type:

```
	make clean
```


# Using `treeducken`

Once `treeducken` has been installed, there are two ways to run the program. The first way entails using command line arguments (see below) to specify the simulation parameters:

* number of simulation replicates (`-r`)
* number of species to simulate (`-nt`)
* species birth rate (`-sbr`)
* species death rate (`-sdr`)
* gene birth rate (`-gbr`)
* gene death rate (`-gdr`)
* gene transfer rate (`-lgtr`)
* individuals to sample per population (`-ipp`)
* effective population size (`-ne`)
* number of loci to simulate on the species tree (`-nl`)
* number of genes to simulate on the locus tree (`-ng`)
* add outgroup with root branch length scaled to fraction (input as fraction) (`-og`)
* outfile prefix (`-o`)
* input settings file (`-i`)


For example you could run:
```
	multiTree -r 10 -nt 100 -sbr 0.5 -sdr 0.4 
```
To simulate 10 replicates of species trees with 100 extant taxa at a species birth rate of 0.5 and a species death rate of 0.2. One may also use the `-i` flag with a settings file (see the example file) that has these parameters (or a combination of the two). 

A bash script, `config_sim.sh`, has been provided that makes the directories where the newick trees are stored and then runs `treeducken` with the `-i simSettings.txt` option. If you would not like the bash script to run `treeducken` simply delete or comment out the line; however, you must run this script to make the correct directory structure to save the tree files. If you do not run this, then the files will not be saved.  

