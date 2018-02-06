# multiTree
multiTree is for the simulation of phylogenies under the three tree model (Rasmussen &amp; Kellis 2012). These three trees represent the species tree, locus (i.e. gene family) tree, and the gene tree.

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

After inputting these parameters the simulation proceeds as follows. First, a species tree is simulated under the constant rate birth-death process until the specified number of species is obtained using the generalized sampling algorithm given in Hartmann 2010. Next, a set of locus tree is forward simulated along the species tree using a birth-death process where transfers are treated as birth events where one loci is transferred to a randomly-selected contemporaneous species, and the other loci stays within the donor lineage. Note that this does simulate along lineages that are bound for extinction. The locus tree simulation process is repeated for the number of loci to simulate on the species tree. For each of those simulated locus trees a gene tree is backward simulated using the multi-locus coalescent. For this multilocus coalescent each tip in the tree starts with the individuals to sample per population. 

# A word of warning

As of February 2018, there has been no unit testing performed. So proceed with caution when simulating...


