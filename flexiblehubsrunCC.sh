#!/bin/bash

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"
subjdir=/camcan/studyforrest/studyforrest/subjlistCC_mo.txt

for chunklen in 15 30 45 60 90 120
do
  for subj in `cat ${subjdir}`
  do
    echo "Requesting chunklen $chunklen for subject $subj"
    sbatch --export=CHUNKLEN=$chunklen,SUBJ=$subj /camcan/FPNflexiblehubs/FPNflexiblehubs/flexiblehubs_1jobCC.sh  
  done
done
