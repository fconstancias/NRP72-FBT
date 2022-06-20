#!/bin/bash

# Parameters:
t=8

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate MetaCHIP

cd /cluster/scratch/haegih/metachip/VAN_CCUG_Donor/

MetaCHIP PI -p Donor -r g -t 8 -i bin -x fa -taxon GTDB.txt
MetaCHIP BP -p Donor -r g -t 8
