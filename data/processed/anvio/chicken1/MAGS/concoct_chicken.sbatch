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