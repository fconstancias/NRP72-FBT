#!/bin/sh

#SBATCH --job-name=anvi-kegg-c
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

anvi-run-kegg-kofams -c CONTIGS.db -T ${OMP_NUM_THREADS}

exit 0
