#!/bin/sh

#SBATCH --job-name=metaph_hum
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

conda activate metaphlan

# JOB STARTS

INPUT=/lustre/constanciasf/NRP72/human/01_QC/
OUT=/lustre/constanciasf/NRP72/human/metaphlan31
mkdir -p ${OUT}


for NAME in `awk '{print $1}' sample_list_metaph.tsv`
do

echo "#### Analyzing " $NAME


ls ${INPUT}${NAME}_R1_host_removed.fastq.gz ${INPUT}${NAME}_R2_host_removed.fastq.gz

echo "#### Starting Default ####"



R1R2s=`ls -m  ${INPUT}${NAME}_R*_host_removed.fastq.gz | sed 's/ //g' | tr -d '\n'`
echo ${R1R2s}

metaphlan ${R1R2s} \
--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 \
--input_type fastq  --bowtie2db /home/constanciasf/db/metaphlan/ \
--force --unknown_estimation \
--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 \
--min_ab 0.000001 \
--legacy-output \
--samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
--add_viruses --tax_lev a --min_cu_len 2000 \
--sample_id ${NAME} --nproc ${OMP_NUM_THREADS} -o ${OUT}/${NAME}_profiled_metagenome.txt

done


# JOB END

merge_metaphlan_tables.py -o ${OUT}/human_merged_profiles.tsv ${OUT}/*_profiled_metagenome.txt

conda activate biobakery3
humann_join_tables  -o ${OUT}/human_merged_profiles.tsv -i ${OUT} --file_name profiled_metagenome.txt

exit 0
