# All agains all binning with SQM pipeline:

_Context:_ 

We generated single assembly MAGs but it is now recognized that all-against-all mapping generate better MAGs.

	SqueezeMeta.pl -m sequential -contiglen 1000 -s /datadrive05/Flo/NRP72/human1/human_SQM_list.tsv -f /datadrive05/Flo/NRP72/	human1/01_QC/ -b 10 -t ${OMP_NUM_THREADS} --nocog --nokegg --nopfam

	SqueezeMeta.pl -m sequential -contiglen 750 -s /datadrive05/Flo/NRP72/chicken1/human_SQM_list.tsv -f /datadrive05/Flo/NRP72/	chicken1/01_QC/ -b 10 -t ${OMP_NUM_THREADS} --nocog --nokegg --nopfam

_Objectives:_ 

- Using the contig files on the SCELSE cluster, loop over the different samples using the all against all strategy.
- Could be a good timing to use another automatic pipeline but well.


	```source ~/.bashrc 
	
	cd /datadrive05/Flo/NRP72/chicken1/ssbins
	
	conda activate /home/constancias/miniconda3/envs/SqueezeMeta260121
	 
	for SAMPLE in `cut  all_against_all.tsv -f1 | sort | uniq`
	do echo "###running " $SAMPLE "###"
	
	ls -lh /datadrive05/Flo/NRP72/chicken1/ssbins/${SAMPLE}/results/01.${SAMPLE}.fasta 
	
	SqueezeMeta.pl -m coassembly -contiglen 1000 --nocog --nokegg --nopfam \
	-p all_against_all_${SAMPLE} \
	-s all_against_all.tsv \
	-extassembly /datadrive05/Flo/NRP72/chicken1/ssbins/${SAMPLE}/results/01.${SAMPLE}.fasta \
	-f /datadrive05/Flo/NRP72/chicken1/01_QC \
	-b 20 -t 10
	
	rm -r all_against_all_${SAMPLE}/data
	rm -r all_against_all_${SAMPLE}/intermediate
	rm -r all_against_all_${SAMPLE}/temp
	done

So it will run in the meantime but no priority for that.