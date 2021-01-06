#!/bin/bash

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"
subjdir=/camcan/studyforrest/studyforrest/subjlistSF.txt

for chunklen in 15 30 45 60 90 120
do
  echo "Requesting chunklen $chunklen for subject $subj"
  sbatch --export=CHUNKLEN=$chunklen,SUBJ=$subj /camcan/studyforrest/studyforrest/flexiblehubs_1jobPoolSubjs.sh 
done
