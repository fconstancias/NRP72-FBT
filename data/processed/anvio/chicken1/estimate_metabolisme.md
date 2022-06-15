	
We upload the ANVIO db to run the `anvi-estimate-metabolism` (anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME):

	DB_files=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/*db
	CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/
	rsync 	-rv -P --checksum ${DB_files} ${CLUSTER_DIR}

		anvi-estimate-metabolism  \
		-p PROFILE.db \
	   -c CONTIGS.db \
	   -C DAS --add-coverage \
	   -O chicken1_estimate_metabo_DAS 
	   
	   
	   	anvi-estimate-metabolism  \
		-p PROFILE.db \
	   -c CONTIGS.db \
	   -C concoct_20 --add-coverage \
	   -O chicken1_estimate_metabo_concoct_20
	   
	   