# activate conda env:

	conda activate anvio-7.1


cd 
	
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
	   
	   
	   anvi-