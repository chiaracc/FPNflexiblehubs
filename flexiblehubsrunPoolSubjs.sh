#!/bin/bash

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"
subjdir=/camcan/studyforrest/studyforrest/subjlistCC.txt

for chunklen in 15 30 45 60 90 120
do
  echo "Requesting chunklen $chunklen"
  sbatch --export=CHUNKLEN=$chunklen /camcan/studyforrest/studyforrest/flexiblehubs_1jobPoolSubjsCC.sh 
done
