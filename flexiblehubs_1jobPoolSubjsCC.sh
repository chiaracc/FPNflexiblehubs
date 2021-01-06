#!/bin/bash
#SBATCH -J FH-CC
#SBATCH --cpus-per-task=8
#SBATCH --output=/camcan/results_rep/camcan/slurm-%j.out
#SBATCH --error=/camcan/results_rep/camcan/slurm-%j.err

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"


echo "In flexiblehubs_1job parameters CHUNKLEN=$CHUNKLEN"
${PYTHON} /camcan/studyforrest/studyforrest/flexibehubsrunCC.py $CHUNKLEN
