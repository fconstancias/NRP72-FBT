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

gtdbtk classify_wf --genome_dir chicken/DAS/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir chicken/DAS/gtdb

gtdbtk classify_wf --genome_dir chicken/concoct_20_refined_2/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir chicken/concoct_20_refined_2/gtdb

gtdbtk classify_wf --genome_dir human/DAS/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir human/DAS/gtdb

gtdbtk classify_wf --genome_dir human/concoct_10_refined_2/ -x fa --cpus ${OMP_NUM_THREADS} --pplacer_cpus ${OMP_NUM_THREADS} --out_dir human/concoct_10_refined_2/gtdb

# JOB END


exit 0
