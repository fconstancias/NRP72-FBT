# Genome-level taxonomic assignments using GTDB:

_Objectives:_ Assign the taxonomy of the MAGs identified using automatic or manual binning.

CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/

##human1:
  
tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/human1-DAS.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/human1-DAS.tar.gz ${CLUSTER_DIR}

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz ${CLUSTER_DIR}

##chicken1:

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz ${CLUSTER_DIR}

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz ${CLUSTER_DIR}

## Transfert:

cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/

ls chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz \
chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz \
human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz \
human/SUMMARY-human/DAS/human1-DAS.tar.gz > list_rsync

rsync -Rdrva --files-from=list_rsync . ${CLUSTER_DIR}

## GTDB-tk


	#!/bin/sh
	
	#SBATCH --job-name=gtdb
	#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
	#SBATCH -N 1 #  Nombre de n<C5><93>uds
	#SBATCH -n 16 # Nombre de t<C3><A2>ches au total
	#SBATCH --mem 100GB # M<C3><A9>moire r<C3><A9>serv<C3><A9>e par n<C5><93>ud
	#SBATCH --mail-user florentinconstancias@gmail.com
	#SBATCH --ntasks-per-node=16
	#SBATCH --ntasks-per-core=1
	#SBATCH --partition=agap_bigmem
	
	export OMP_NUM_THREADS=16
	
	module load python/Anaconda/3-5.1.0
	eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980
	
	conda activate /home/constanciasf/.conda/envs/gtdbtk
	# JOB STARTS
	
	cd /home/constanciasf/scratch/NRP72/gtdb2.1/genomes
	
	gtdbtk classify_wf --genome_dir chicken/DAS/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir chicken/DAS/gtdb --write_single_copy_genes chicken/DAS/gtdb/single_copy_genes.fa
	
	gtdbtk classify_wf --genome_dir chicken/concoct_20_refined_2/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir chicken/concoct_20_refined_2/gtdb --write_single_copy_genes chicken/concoct_20_refined_2/gtdb/single_copy_genes.fa
	
	gtdbtk classify_wf --genome_dir human/DAS/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir human/DAS/gtdb --write_single_copy_genes human/DAS/gtdb/single_copy_genes.fa
	
	gtdbtk classify_wf --genome_dir human/concoct_10_refined_2/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir human/concoct_10_refined_2/gtdb --write_single_copy_genes human/concoct_10_refined_2/gtdb/single_copy_genes.fa
	
	# JOB END
	
	
	exit 0