Install:
https://github.com/gem-pasteur/Integron_Finder

doc:
https://integronfinder.readthedocs.io/en/v2.0.2/user_guide/tutorial.html

cd /datadrive05/Flo/tools/PathoFact_orf

conda activate integron_finder
integron_finder chicken_l1000.fna --cpu 20 --calin-threshold 1 --func-annot --linear --pdf --gbk --path_func_annot bank_hmm #--promoter-attI

able to find :NCBIfam-AMRFinder.hmm ??

find . -name "NCBIfam-AMRFinder.hmm" -print 

vim bank_hmm
/home/constancias/miniconda3/envs/atlas-metagenomes/envs/integron_finder/lib/python3.9/site-packages/integron_finder/data/Functional_annotation/NCBIfam-AMRFinder.hmm
/home/constancias/miniconda3/envs/atlas-metagenomes/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/Pfam/Pfam-A.hmm
/home/constancias/miniconda3/envs/atlas-metagenomes/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/KEGG/HMMs/Kofam.hmm

integron_finder --cpu 20 --calin-threshold 1 --func-annot --linear --pdf --gbk \
--path-func-annot bank_hmm --outdir Results_Integron_Finder_chicken_l1000_2 chicken_l1000.fna