# Preparation:

## First, we rename the genomes based on samples:

	cd /datadrive05/Flo/NRP72/chicken1/ssbins/
	for d in *; do 
	cd ${d}/results/DAS/*ASTool_bins/
	for f in *fa
	#do echo ${d}; ls "${f}"
	do mv -- "${f}" "${d}_${f}"
	done
	cd /datadrive05/Flo/NRP72/chicken1/ssbins/
	done

## We also store those in a tar.gz archive

mkdir ssbins_genomes_chicken1
scp */results/DAS/*_DASTool_bins/*fa ssbins_genomes_chicken1
cd ssbins_genomes_chicken1
for f in *fa
	do mv -- "${f}" "chicken1_${f}"
	done
cd ..
tar -czvf ssbins_genomes_chicken1.tar.gz ssbins_genomes_chicken1

# dRep :
source ~/.bashrc 
conda activate  /home/constancias/miniconda3/envs/drep
	 ...::: dRep v3.2.0 :::...

## only SS bins:

### 0.95
dRep dereplicate chicken1_drep_SS_095 \
-g ssbins_genomes_chicken1/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.95 \


### 0.98
dRep dereplicate chicken1_drep_SS_098 \
-g ssbins_genomes_chicken1/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.98 \

## gtdb
conda activate gtdbtk

gtdbtk classify_wf --genome_dir drep_all/dereplicated_genomes/ -x fa --cpus 10 --pplacer_cpus 10 --out_dir drep_all/dereplicated_genomes/gtdb


## SS bins + coassembly bins:

