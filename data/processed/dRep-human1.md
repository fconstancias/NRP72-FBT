# Preparation:

First, we rename the genomes based on samples:

	cd /datadrive05/Flo/NRP72/human1/ssbins/
	for d in *; do 
	cd ${d}/results/DAS/*ASTool_bins/
	for f in *fa
	#do echo ${d}; ls "${f}"
	do mv -- "${f}" "${d}_${f}"
	done
	cd /datadrive05/Flo/NRP72/human1/ssbins/
	done

## We also store those in a tar.gz archive

mkdir ssbins_genomes_human1
scp */results/DAS/*_DASTool_bins/*fa ssbins_genomes_human1
cd ssbins_genomes_human1
for f in *fa
	do mv -- "${f}" "human1_${f}"
	done
cd ..
tar -czvf ssbins_genomes_human1.tar.gz ssbins_genomes_human1

# dRep :
source ~/.bashrc 
conda activate  /home/constancias/miniconda3/envs/drep
 	...::: dRep v3.2.0 :::...
## only SS bins:

### 0.95
dRep dereplicate human1_drep_SS_095 \
-g ssbins_genomes_human1/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.95 \


### 0.98
dRep dereplicate human1_drep_SS_098 \
-g ssbins_genomes_human1/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.98 \



## coassembly bins:

cd /datadrive05/Flo/NRP72/dRep_SS_coASS/human1

### 0.95

dRep dereplicate human1_drep_DAS_concoct_10_refined_2_095 \
-g DAS/*fa concoct_10_refined_2/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.95

### 0.98

dRep dereplicate human1_drep_DAS_concoct_10_refined_2_098 \
-g DAS/*fa concoct_10_refined_2/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.98

## coassembly + ss bins:

### DAS

### 0.95

dRep dereplicate human1_drep_DAS_SS_095 \
-g DAS/*fa SS/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.95

### 0.98

dRep dereplicate human1_drep_DAS_SS_098 \
-g DAS/*fa SS/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.98

### manual

### 0.95

dRep dereplicate human1_drep_concoct_10_refined_2_SS_095 \
-g concoct_10_refined_2/*fa SS/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.95

### 0.98

dRep dereplicate human1_drep_concoct_10_refined_2_SS_098 \
-g concoct_10_refined_2/*fa SS/*fa \
-p 10 -comp 70 -con 20 \
--checkM_method lineage_wf --S_algorithm ANImf \
-sa 0.98

