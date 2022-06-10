	
We upload the ANVIO db to run the `anvi-estimate-metabolism` (anvi-estimate-metabolism -c CONTIGS.db -p PROFILE.db -C COLLECTION_NAME):

	DB_files=/Users/fconstan/Documents/GitHub/NRP72-FBT/data/processed/anvio/chicken/*db
	CLUSTER_DIR=constanciasf@muse-login.meso.umontpellier.fr:/lustre/constanciasf/NRP72/chicken/ANVIO_WDIR/
	rsync 	-rv -P --checksum ${DB_files} ${CLUSTER_DIR}
