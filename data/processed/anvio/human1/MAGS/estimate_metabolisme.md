
We upload the ANVIO db to run the `anvi-run-kegg-kofams` and `anvi-estimate-metabolism`:

	DB_files=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/human/*db
	CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/human/ANVIO_WDIR/
	
	rsync 	-rv -P --checksum ${DB_files} ${CLUSTER_DIR}

	anvi-run-kegg-kofams -c CONTIGS.db -T ${OMP_NUM_THREADS}

	
		anvi-estimate-metabolism  \
		-p PROFILE.db \
	   -c CONTIGS.db \
	   -C DAS --add-coverage \
	   -O human_estimate_metabo_DAS 
	   
	   
	   	anvi-estimate-metabolism  \
		-p PROFILE.db \
	   -c CONTIGS.db \
	   -C concoct_10 --add-coverage \
	   -O human_estimate_metabo_concoct_10