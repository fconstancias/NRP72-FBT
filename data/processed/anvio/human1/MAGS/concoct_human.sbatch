#!/bin/sh

#SBATCH --job-name=human_concoct
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

CONT_DB=/lustre/constanciasf/NRP72/human/ANVIO_WDIR/CONTIGS.db
PROF=/lustre/constanciasf/NRP72/human/ANVIO_WDIR/PROFILE.db

for nclust in 20 30 50 70
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

exit 0