#!/bin/sh

#SBATCH --job-name=plasx
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -N 1 #  Nombre de nœuds
#SBATCH -n 10 # Nombre de tâches au total
#SBATCH --mem 100GB # Mémoire réservée par nœud
#SBATCH --mail-user h.haegi@gmail.com
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --partition=agap_normal

export OMP_NUM_THREADS=10

module load python/Anaconda/3-5.1.0
eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980


PREFIX='chicken'
THREADS=10

cd /home/constanciasf/scratch/final_nrp72/raw/

# Activate anvio
conda activate anvio-7.1

# Create anvio contigs db and export gene calls
anvi-gen-contigs-database -L 0 -T $THREADS --project-name $PREFIX -f $PREFIX.fa -o $PREFIX.db

anvi-export-gene-calls --gene-caller prodigal -c $PREFIX.db -o $PREFIX-gene-calls.txt


#Download COGs and Pfams db
#anvi-setup-ncbi-cogs --cog-version COG14 --cog-data-dir COG_2014 -T $THREADS
#anvi-setup-pfams --pfam-version 32.0 --pfam-data-dir Pfam_v32

#Annotate COGs and Pfams and export to text files
anvi-run-ncbi-cogs -T $THREADS --cog-version COG14 --cog-data-dir COG_2014 -c $PREFIX.db
anvi-run-pfams -T $THREADS --pfam-data-dir Pfam_v32 -c $PREFIX.db
anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c $PREFIX.db -o $PREFIX-cogs-and-pfams.txt

#Annotate de novo gene families
conda deactivate
conda activate plasx
#plasx setup \
#    --de-novo-families 'https://zenodo.org/record/5819401/files/PlasX_mmseqs_profiles.tar.gz' \
#    --coefficients 'https://zenodo.org/record/5819401/files/PlasX_coefficients_and_gene_enrichments.txt.gz'

#FYI: the command “plasx setup” needs to be run the first time, in the script I send you it’s marked as comment because I had some trouble running the script and run it therefore several times, but since the setup was already done, I did not need to run that again (it also takes some while to run)
    
#Annotate genes using the database of de novo gene families
plasx search_de_novo_families \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-de-novo-families.txt \
    --threads $THREADS \
    --splits 32 \
    --overwrite
    
    
#Use PlasX to classify contigs as plasmid or non-plasmid sequences
plasx predict \
    -a $PREFIX-cogs-and-pfams.txt $PREFIX-de-novo-families.txt \
    -g $PREFIX-gene-calls.txt \
    -o $PREFIX-scores.txt \
    --overwrite
