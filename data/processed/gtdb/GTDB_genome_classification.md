# Genome-level taxonomic assignments using GTDB:

_Objectives:_ Assign the taxonomy of the MAGs identified using automatic or manual binning.

CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/

##human1:
  
tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/human1-DAS.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/DAS/human1-DAS.tar.gz ${CLUSTER_DIR}

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz ${CLUSTER_DIR}

##chicken1:

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz ${CLUSTER_DIR}

tar -czf /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz \
/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/bin_by_bin/*/*-contigs.fa

rsync 	-rv -P /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz ${CLUSTER_DIR}

## Transfert:

cd /Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/

ls chicken/SUMMARY-chicken/concoct_20_refined_2/chicken1-concoct_20_refined_2.tar.gz \
chicken/SUMMARY-chicken/DAS/chicken1-DAS.tar.gz \
human/SUMMARY-human/concoct_10_refined_2/human1-concoct_10_refined_2.tar.gz \
human/SUMMARY-human/DAS/human1-DAS.tar.gz > list_rsync

rsync -Rdrva --files-from=list_rsync . ${CLUSTER_DIR}

## GTDB-tk

