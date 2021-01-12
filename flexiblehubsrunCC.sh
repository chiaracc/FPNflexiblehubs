#!/bin/bash

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"
subjdir=/camcan/FPNflexiblehubs/FPNflexiblehubs/subjlistCC.txt

for chunklen in 15 30 45 60 90 120
do
  echo "Requesting chunklen $chunklen"
  sbatch --export=CHUNKLEN=$chunklen /camcan/FPNflexiblehubs/FPNflexiblehubs/flexiblehubs_1jobCC.sh 
done
