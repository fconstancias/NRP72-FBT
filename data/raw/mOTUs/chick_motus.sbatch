#!/bin/sh

#SBATCH --job-name=chick_motus
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -N 1 #  Nombre de noeuds
#SBATCH -n 10 # Nombre de tacches au total
#SBATCH --mem 8GB # Memoire reservee par noeud
#SBATCH --mail-user florentinconstancias@gmail.com
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --partition=agap_normal

export OMP_NUM_THREADS=10

module load python/Anaconda/3-5.1.0 
conda --version

eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980

conda activate mOTUs_v2

READS_DNA=/lustre/constanciasf/NRP72/chicken/01_QC/
MOTUS=/lustre/constanciasf/NRP72/chicken/mOTUs_v2/
MOTUS_BAM=/lustre/constanciasf/NRP72/chicken/mOTUs_v2/BAM/
MOTUS_RC=/lustre/constanciasf/NRP72/chicken/mOTUs_v2/RC/
MOTUS_PROF=/lustre/constanciasf/NRP72/chicken/mOTUs_v2/PROFILE/

mkdir -p ${MOTUS_PROF}; mkdir -p ${MOTUS_RC}; mkdir -p ${MOTUS_BAM}

for sample in `awk '{print $1}' sample_list_metaph.tsv`
do

ls -lh ${READS_DNA}/${sample}_R1_diet_removed.fastq.gz ${READS_DNA}/${sample}_R2_diet_removed.fastq.gz

/home/constanciasf/tools/mOTUs_v2/motus profile \
    -f ${READS_DNA}/${sample}_R1_diet_removed.fastq.gz \
    -r ${READS_DNA}/${sample}_R2_diet_removed.fastq.gz \
    -t ${OMP_NUM_THREADS} \
    -n ${sample} \
    -I ${MOTUS_BAM}${sample}.bam \
    -M ${MOTUS_RC}${sample}_mgc_read_count.txt \
    -q -u \
    -g 3 \
    -l 75 \
    -o ${MOTUS_PROF}${sample}_profile.txt

done

# If you want to run mOTUs with different parameters, you can save the intermediate bam file and provide it as input afterwards. Example:
# https://github.com/motu-tool/mOTUs_v2/wiki/Profile-one-sample#save-profile-to-file

#motu profile -s test1.fq -I test1.bam -l 40 -o result_l40.motus
#motu profile -i test1.bam -l 70 -o result_l70.motus
#motu profile -i test1.bam -l 100 -o result_l100.motus

/home/constanciasf/tools/mOTUs_v2/motus merge -d ${MOTUS_PROF} -o ${MOTUS_PROF}merged_profile.tsv
# /homedir/constancias/tools/mOTUs_v2/motus merge -d ${MOTUS_PROF} -a all -o ${MOTUS_PROF}merged_with_all_profile.tsv

exit 0


