#!/bin/bash
touch speciestrees.tre
touch locustrees.tre
touch genetrees.tre
for((i = 0; i < 10; ++i)) 
do
	cat ./sim_files/speciestree_$i/*.sp.tre >> ./speciestrees.tre
	echo "" >> ./speciestrees.tre
	for((j = 0; j < 10; ++j))
	do
		cat ./sim_files/speciestree_$i/*_$i_$j.loc.tre >> ./locustrees_$i.tre
		echo "" >> ./locustrees.tre
		cat ./sim_files/speciestree_$i/*_$i_$j.gen.tre >> ./genetrees_$i.tre
		echo "" >> ./genetrees.tre
	done
done

