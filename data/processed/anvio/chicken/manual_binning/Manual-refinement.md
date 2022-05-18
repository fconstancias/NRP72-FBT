
# Manual binning of the chicken metagenomes with Anvi'o

SQM pipeline generated lots of data including MAGs binned based on automatic approach. We are going to adopt here the manual binning approach starting from CONCOCT clusters.

We first upload the .db files to a cluster to perform the clustering using CONCOCT.

	DB_files=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/*db
	CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/
	SLOTS=6
	
	rsync 	-rv -P --checksum ${DB_files} ${CLUSTER_DIR}
	
Then we use anvio-7.1 to perform clustering with CONCOCT

	#!/bin/sh
	
	#SBATCH --job-name=chick_concoct
	#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
	#SBATCH -N 1 #  Nombre de nœuds
	#SBATCH -n 10 # Nombre de tâches au total
	#SBATCH --mem 100GB # Mémoire réservée par nœud
	#SBATCH --mail-user florentinconstancias@gmail.com
	#SBATCH --ntasks-per-node=10
	#SBATCH --ntasks-per-core=1
	#SBATCH --partition=agap_normal
	
	export OMP_NUM_THREADS=10
	
	module load python/Anaconda/3-5.1.0
	eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980
	
	
	# JOB STARTS
	
	
	conda activate anvio-7.1
	
	CONT_DB=/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/CONTIGS.db
	PROF=/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/PROFILE.db
	
	for nclust in 20
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
	
	anvi-split -p ${PROF} \
	-c ${CONT_DB} \
	-C concoct_20 \
	-o split_concoct_20 --skip-variability-tables --skip-hierarchical-clustering
	
	anvi-cluster-contigs \
	-p split_concoct_20/Bin_5/PROFILE.db \
	-c split_concoct_20/Bin_5/CONTIGS.db \
	-T ${OMP_NUM_THREADS} \
	--driver concoct \
	--clusters 4 \
	-C concoct_20_Bin_5_3 --just-do-it
	
	anvi-summarize \
	-p split_concoct_20/Bin_5/PROFILE.db \
	-c split_concoct_20/Bin_5/CONTIGS.db \
	-C concoct_20_Bin_5_3 \
	-o SUMMARY/concoct_20_Bin_5_3
	
	anvi-export-collection -C concoct_20_Bin_5_3 \
	-p split_concoct_20/Bin_5/PROFILE.db \
	 -O concoct_concoct_20_Bin_5_3_collection.txt
	
	exit 0


We download the generated collections.

	WD_DIR_LOCAL=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/manual_binning/
	CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/
	
	rsync 	-rv -P --checksum ${CLUSTER_DIR}*_collection.txt.txt ${WD_DIR_LOCAL}
	rsync -rv -P --relative --checksum ${CLUSTER_DIR}SUMMARY/*/bins_summary.txt ${WD_DIR_LOCAL}

We will try to start with the `concoct_20` collection since eventough bin_0 = 28836 bins. We will see :

	conda activate anvio-7.1
	DB_dir=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/
	cd ${DB_dir}
	PROF=PROFILE.db
	CONT_DB=CONTIGS.db
	
	anvi-summarize -c ${CONT_DB} -p ${PROF} --list-collections
	
	anvi-import-collection  -c ${CONT_DB} -p ${PROF} -C concoct_20 manual_binning/concoct_20_collection.txt.txt 
	
We are going to refine this freshly imported concoct_20 collection. First we bin roughly we will go for a second refinement later.

	for bin in `cut -f2  manual_binning/concoct_20_collection.txt.txt | uniq`
	do
	echo "### Refining bin" ${bin} "Start###"
	
	        anvi-refine -c ${CONT_DB} \
	                    -p ${PROF}  \
	                    -C concoct_20 \
	                    -b ${bin}
	done

	anvi-export-collection -C concoct_20 \
	-p ${PROF} \
	-O manual_binning/concoct_20_collection_refined
	
	anvi-import-collection  -c ${CONT_DB} -p ${PROF} -C concoct_20_raw manual_binning/concoct_20_collection.txt.txt 

	anvi-summarize -c ${CONT_DB} -p ${PROF} -C concoct_20 -o manual_binning/SUMMARY/concoct_20_refined

Start diff coverage and seq composition in detection - bin - adapt with only diff coverage and some times then only seq composition - when coverage in only one sample.