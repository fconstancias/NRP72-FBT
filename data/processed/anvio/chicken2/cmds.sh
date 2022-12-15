source ~/.bashrc 


#
cd /datadrive05/Flo/NRP72/01_RAW/2/

ls *.fq.gz | cut -d_ -f1,2 | sort | uniq | while read SAMPLE; do cat ${SAMPLE}*_1.fq.gz > ${SAMPLE}_R1_.fastq.gz && cat ${SAMPLE}*_2.fq.gz > ${SAMPLE}_R2_.fastq.gz ; done

#
conda activate kneaddata

cd /datadrive05/Flo/db/Bowtie2Index/

kneaddata_database --download human_genome bowtie2 /datadrive05/Flo/db/Bowtie2Index/human_genome_kneadd/

wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz .

bowtie2-build --threads 6 --seed 123 GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz  CF_016699485.2_bGalGal1

wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz  .

bowtie2-build --threads 6 --seed 123 GCA_000786075.2_hs38d1_genomic.fna.gz GCA_000786075.2_hs38d1_genomic

bowtie2-build --threads 6 --seed 123 genome.fa PhiX_kneaddata

#
cd /datadrive05/Flo/NRP72/01_RAW/2/

ls *gz | cut -d_ -f1,2  | sort | uniq > ../../samples_list_2.tsv


cd /datadrive05/Flo/NRP72/

INPUT=/datadrive05/Flo/NRP72/01_RAW/2/
OUTPUT=/datadrive05/Flo/NRP72/kneaddata_out/

for SAMPLE in D_C; do

echo "###running Sample" ${SAMPLE}  "###"

kneaddata -i ${INPUT}${SAMPLE}_R1_.fastq.gz -i ${INPUT}${SAMPLE}_R2_.fastq.gz \
-db /datadrive05/Flo/db/Bowtie2Index/CF_016699485.2_bGalGal1 \
--run-trf --run-trim-repetitive --decontaminate-pairs strict \
--run-fastqc-start --run-fastqc-end --output-prefix ${SAMPLE} --remove-intermediate-output \
-t 4 \
--trimmomatic /home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/ \
--trimmomatic-options="ILLUMINACLIP:/home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:3  MINLEN:60" \
-o ${OUTPUT} 

done

for SAMPLE in D_H; do

echo "###running Sample" ${SAMPLE}  "###"

kneaddata -i ${INPUT}${SAMPLE}_R1_.fastq.gz -i ${INPUT}${SAMPLE}_R2_.fastq.gz \
-db /datadrive05/Flo/db/Bowtie2Index/GCA_000786075.2_hs38d1_genomic \
--run-trf --run-trim-repetitive --decontaminate-pairs strict \
--run-fastqc-start --run-fastqc-end --output-prefix ${SAMPLE} --remove-intermediate-output \
-t 4 \
--trimmomatic /home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/ \
--trimmomatic-options="ILLUMINACLIP:/home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:3  MINLEN:60" \
-o ${OUTPUT} 

done

ls ${OUTPUT}/*paired_{1,2}.fastq | parallel -j 4 -t gzip


cd /datadrive05/Flo/NRP72/

INPUT=/datadrive05/Flo/NRP72/01_RAW/2/
SAMPLELIST=/datadrive05/Flo/NRP72/samples_list_2_noDonors.tsv
OUTPUT=/datadrive05/Flo/NRP72/kneaddata_out_nodonors/

#for SAMPLE in `awk '{print $1}' ${SAMPLELIST}`; do
for SAMPLE in C7_8 ; do 

####-db /datadrive05/Flo/db/Bowtie2Index/PhiX/PhiX_kneaddata
echo "###running Sample" ${SAMPLE}  "###"

kneaddata -i ${INPUT}${SAMPLE}_R1_.fastq.gz -i ${INPUT}${SAMPLE}_R2_.fastq.gz \
-db /datadrive05/Flo/db/Bowtie2Index/PhiX/PhiX_kneaddata \
--run-trf --run-trim-repetitive --decontaminate-pairs strict \
--run-fastqc-start --run-fastqc-end --output-prefix ${SAMPLE} --remove-intermediate-output \
-t 4 \
--trimmomatic /home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/ \
--trimmomatic-options="ILLUMINACLIP:/home/constancias/miniconda3/envs/atlas-metagenomes/envs/kneaddata/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:3  MINLEN:60" \
-o ${OUTPUT} 


ls ${OUTPUT}/${SAMPLE}*paired_{1,2}.fastq | parallel -j 2 -t gzip

done

# coassembly all cleaned reads:

conda activate  /home/constancias/miniconda3/envs/atlas-metagenomes/envs/anvio-7.1

NSLOTS=10
IN=/datadrive05/Flo/NRP72/kneaddata_out/chicken2
OUT=/datadrive05/Flo/NRP72/02_ASSEMBLY/

mkdir -p ${OUT}
MIN_CONTIG_SIZE=1000

R1s=`ls -m  ${IN}/*paired_1.fastq.*gz | sed 's/ //g' | tr -d '\n'`
R2s=`ls -m  ${IN}/*paired_2.fastq.*gz | sed 's/ //g' | tr -d '\n'`


megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.20 -o ${OUT}/chicken2 -t $NSLOTS
# did not work.

#
conda activate /home/constancias/miniconda3/envs/SqueezeMeta160321

cd /datadrive05/Flo/NRP72/SQM

SqueezeMeta.pl -m coassembly -contiglen 750 --cleaning --nobins \
-s /datadrive05/Flo/NRP72/sample_list_chicken2.tsv \
-f /datadrive05/Flo/NRP72/01_RAW/2/chicken2/ \
-b 12 -t 12 -p chicken2

# parallel -j 8 'samtools view -bS {} -@ 4 > {}.bam' ::: *.sam 
# rename  's/.bam/-RAW.bam/' *bam
# mv *bam ../../results/sqm2anvio/bam
# rename  's/chicken2.//' *bam

source /home/constancias/miniconda3/etc/profile.d/conda.sh

conda activate anvio-7.1

python3 /home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py \
-p chicken2 -o chicken2/results/anvio --num-threads 10 --min-contig-length 1000 --run-HMMS  --run-scg-taxonomy --profile-SCVs 

#humann2                                                                                                                                                                      
python3 /home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py \
-p human2 -o human2/results/anvio --num-threads 10 --min-contig-length 1000 --run-HMMS  --force-overwrite  --profile-SCVs 



anvi-run-hmms -c chicken2/results/anvio/CONTIGS.db -T 6 --also-scan-trnas
anvi-run-scg-taxonomy -c chicken2/results/anvio/CONTIGS.db -T 4
anvi-run-kegg-kofams -c chicken2/results/anvio/CONTIGS.db -T 6

# metabin
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

# from there we check the SUMMARY/concoct_${nclust}/bins_summary.tsv
# 20 clusters is the way to go bins <25000 splits BUT Bin_4 (> 60 000 contigs) has to be further clustered.
# still there with 50 metabins...

anvi-split  -p ${PROF} \
-c ${CONT_DB} \
-C concoct_20 \
-b Bin_4 \
-o split_concoct_20 --skip-variability-tables --skip-hierarchical-clustering

# redefine the .db for the metabin impossible to bin:
PROF=split_concoct_20/Bin_4/PROFILE.db
CONT_DB=split_concoct_20/Bin_4/CONTIGS.db

# run concoct on this particular bin:

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



anvi-run-hmms -c human2/results/anvio/CONTIGS.db -T 6 --also-scan-trnas --just-do-it
anvi-run-scg-taxonomy -c human2/results/anvio/CONTIGS.db -T 4
anvi-run-kegg-kofams -c human2/results/anvio/CONTIGS.db -T 6



✓ anvi-run-hmms took 0:31:18.512131
diamond v2.0.8.146 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

diamond v2.0.8.146 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org

Traceback (most recent call last):
  File "/home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py", line 433, in <module>
    main(parse_arguments())
  File "/home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py", line 412, in main
    create_contigsDB(args.project, contigs, genes, functions, tax_genes, tax_contigs, args.outputDir, args.run_HMMS, args.num_threads, args.run_scg_taxonomy, version, logfile)
  File "/home/constancias/miniconda3/envs/SqueezeMeta160321/bin/anvi-load-sqm.py", line 228, in create_contigsDB
    if diamond_version_path != diamond_version_makedb:
UnboundLocalError: local variable 'diamond_version_path' referenced before assignment
(anvio-7.1) constancias@scelse:/datadrive05/Flo/NRP72/SQM$ 

#
conda activate /home/constancias/miniconda3/envs/SqueezeMeta160321

cd /datadrive05/Flo/NRP72/SQM

SqueezeMeta.pl -m coassembly -contiglen 750 --cleaning --nobins \
-s /datadrive05/Flo/NRP72/sample_list_human2.tsv \
-f /datadrive05/Flo/NRP72/01_RAW/2/human2/ \
-b 12 -t 12 -p human2


# mobileOG-db:

diamond makedb --in mobileOG-db_beatrix-1.6.All.faa  -d mobileOG-db_beatrix-1.6.All.dmnd --threads 6
#chmod +x mobileOGs-pl-kyanite.sh
DIAMOND=/datadrive05/Flo/db/mobileOG-db_workdir/mobileOG-db_beatrix-1.6.All.dmnd
KVALUE=15
ESCORE=1e-20
QUERYSCORE=90
PIDENTVALUE=90
METADATA=/datadrive05/Flo/db/mobileOG-db_workdir/mobileOG-db-beatrix-1.6-All.csv
#in_faa=/datadrive05/Flo/NRP72/SQM/chicken2/results/03.chicken2
#in_faa=/datadrive05/Flo/NRP72/SQM/human2/results/03.human2
#in_faa=/datadrive05/Flo/db/mobileOG-db_workdir/03.chicken
#in_faa=/datadrive05/Flo/db/mobileOG-db_workdir/03.SqueezeHuman

#out=/datadrive05/Flo/db/mobileOG-db_workdir/03.chicken2
#out=/datadrive05/Flo/db/mobileOG-db_workdir/03.human2
#out=/datadrive05/Flo/db/mobileOG-db_workdir/03.chicken1
#out=/datadrive05/Flo/db/mobileOG-db_workdir/03.SqueezeHuman

diamond blastp -q ${in_faa}.faa --db ${DIAMOND} --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend -k $KVALUE -o ${out}.tsv -e $ESCORE --query-cover $QUERYSCORE --id $PIDENTVALUE --threads 4
python mobileOGs-pl-kyanite.py --o ${out} --i ${out}.tsv -m ${METADATA}

#geNomad

conda activate genomad

#genomad download-database .
#in=/datadrive05/Flo/NRP72/SQM/human2/results/01.human2.fasta
#out=/datadrive05/Flo/db/genomad_test/03.human2

#in=/datadrive05/Flo/db/genomad_test/01.SqueezeHuman.fasta
#out=/datadrive05/Flo/db/genomad_test/human1

#in=/datadrive05/Flo/NRP72/SQM/chicken2/results/01.chicken2.fasta
#out=/datadrive05/Flo/db/genomad_test/chicken2

#in=/datadrive05/Flo/db/genomad_test/01.chicken.fasta 
#out=/datadrive05/Flo/db/genomad_test/chicken1


cd /datadrive05/Flo/db/genomad_test
genomad end-to-end --min-score 0.7 --splits 2 ${in} ${out} /datadrive05/Flo/db/genomad_db -t 4



