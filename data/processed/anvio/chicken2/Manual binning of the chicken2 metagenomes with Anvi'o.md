
# Manual binning of the chicken metagenomes with Anvi'o

SQM pipeline generated lots of data including MAGs binned based on automatic approach. We are going to adopt here the manual binning approach starting from CONCOCT clusters.

We first upload the .db files to a cluster to perform the clustering using CONCOCT.

	source /home/constancias/miniconda3/etc/profile.d/conda.sh
	
	conda activate anvio-7.1
	
	python3 /home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py \
	-p chicken2 -o chicken2/results/anvio --num-threads 10 --min-contig-length 1000 --run-HMMS  --run-scg-taxonomy --profile-SCVs 

	
	anvi-run-hmms -c chicken2/results/anvio/CONTIGS.db -T 6 --also-scan-trnas
	anvi-run-scg-taxonomy -c chicken2/results/anvio/CONTIGS.db -T 4
	anvi-run-kegg-kofams -c chicken2/results/anvio/CONTIGS.db -T 6


Then we use anvio-7.1 to perform clustering with CONCOCT

	
#metabin
	
	cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken2/
	PROF=PROFILE.db
	CONT_DB=CONTIGS.db
	OMP_NUM_THREADS=6
	
	for nclust in 7 10 15 25 40 50
	do
	
	anvi-cluster-contigs \
	-p ${PROF} \
	-c ${CONT_DB} \
	-T ${OMP_NUM_THREADS} \
	--driver concoct \
	--clusters ${nclust} \
	-C concoct_${nclust} --just-do-it
	
	anvi-summarize \
	-p ${PROF} \
	-c ${CONT_DB} \
	-C concoct_${nclust} \
	-o SUMMARY/concoct_${nclust}
	
	anvi-export-collection -C concoct_${nclust} \
	-p ${PROF} \
	-O concoct_${nclust}_collection.txt
	
	done

from there we check the SUMMARY/concoct_${nclust}/bins_summary.tsv
20 clusters is the way to go bins <25000 splits BUT Bin_4 (> 60 000 contigs) has to be further clustered.
still there with 50 metabins...



	anvi-split  -p ${PROF} \
	-c ${CONT_DB} \
	-C concoct_20 \
	-b Bin_4 \
	-o split_concoct_20 --skip-variability-tables --skip-hierarchical-clustering
	

# redefine the .db for the metabin impossible to bin:
PROF=split_concoct_20/Bin_4/PROFILE.db
CONT_DB=split_concoct_20/Bin_4/CONTIGS.db


run concoct on this particular bin:

	for nclust in 2 3 5 10
	do
	anvi-cluster-contigs \
	    -p ${PROF}  \
	    -c ${CONT_DB}  \
	    -T ${OMP_NUM_THREADS} \
	    --driver concoct \
	    --clusters ${nclust} \
	    -C concoct_20_Bin_4_${nclust} --just-do-it
	
	anvi-summarize \
	    -p ${PROF}  \
	    -c ${CONT_DB}  \
	    -C concoct_20_Bin_4_${nclust} \
	    -o SUMMARY/concoct_20_Bin_4_${nclust}
	
	anvi-export-collection -C concoct_20_Bin_4_${nclust} \
	    -p ${PROF}  \
	    -O concoct_20_Bin_4_${nclust}_collection

	done

Again, we explore the SUMMARY/concoct_20_Bin_4_${nclust}/bins_summary.txt. 5 clusters for this bin is the way to go!

We need to combine the concoct_20_collection.txt.txt and the concoct_20_Bin_4_5_collection.txt

So we need to remove Bin_4 from concoct_20_collection.txt.txt and concatenate with concoct_20_Bin_4_5_collection.txt.

	cut -f2  concoct_20_collection.txt.txt | uniq
	grep -v "\<Bin_4\>" concoct_20_collection.txt.txt | cut -f2  | uniq
	sed 's/Bin_/Bin_4_/g' concoct_20_Bin_4_5_collection.txt | cut -f2 | uniq
	
	grep -v "\<Bin_4\>" concoct_20_collection.txt.txt  > concoct_20_collection_nobin4.tsv
	sed 's/Bin_/Bin_4_/g' concoct_20_Bin_4_5_collection.txt > concoct_20_Bin_4_5_collection_renamed.tsv
	cat concoct_20_collection_nobin4.tsv concoct_20_Bin_4_5_collection_renamed.tsv > concoct_20_20_Bin_4_concoct_5_merged.tsv
	
	
An we then import that in the chicken2 databases.


	cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken2/
	PROF=PROFILE.db
	CONT_DB=CONTIGS.db
	anvi-import-collection  -c ${CONT_DB} -p ${PROF} -C concoct_20_Bin_4_concoct_5 concoct_20_20_Bin_4_concoct_5_merged.tsv



We are going to refine this freshly imported concoct_20_Bin_4_concoct_5 collection. First we bin roughly we will go for a second refinement later.

	for bin in `cut -f2 concoct_20_20_Bin_4_concoct_5_merged.tsv | uniq`
	do
	echo "### Refining bin" ${bin} "Start###"
	
	        anvi-refine -c ${CONT_DB} \
	                    -p ${PROF}  \
	                    -C concoct_20_Bin_4_concoct_5 \
	                    -b ${bin}
	done



	anvi-export-collection -C concoct_20 \
	-p ${PROF} \
	-O manual_binning/concoct_20_collection_refined
	
	anvi-import-collection  -c ${CONT_DB} -p ${PROF} -C concoct_20_raw manual_binning/concoct_20_collection.txt.txt 

	anvi-summarize -c ${CONT_DB} -p ${PROF} -C concoct_20 -o manual_binning/SUMMARY/concoct_20_refined

Start diff coverage and seq composition in detection - bin - adapt with only diff coverage and some times then only seq composition - when coverage in only one sample.

We now go for a second run of binning: refinement!

	for bin in `cut -f2  manual_binning/concoct_20_collection_refined.txt| uniq`
	do
	echo "### Refining bin" ${bin} "Start###"
	
	        anvi-refine -c ${CONT_DB} \
	                    -p ${PROF}  \
	                    -C concoct_20 \
	                    -b ${bin}
	done
	
	anvi-export-collection -C concoct_20 \
	-p ${PROF} \
	-O manual_binning/concoct_20_collection_refined_2

	anvi-summarize -c ${CONT_DB} -p ${PROF} -C concoct_20 -o SUMMARY/concoct_20_refined_2



Export the DAS collection summary including genomes as well as refined2:

	conda activate anvio-7.1
	DB_dir=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/
	cd ${DB_dir}
	PROF=PROFILE.db
	CONT_DB=CONTIGS.db
	
	anvi-summarize -c ${CONT_DB} -p ${PROF} --list-collections

	anvi-summarize -c ${CONT_DB} -p ${PROF} -C DAS -o SUMMARY/DAS/
	
	anvi-summarize -c ${CONT_DB} -p ${PROF} -C concoct_20 -o SUMMARY/concoct_20_refined_2
	
