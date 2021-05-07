#!/bin/bash
#SBATCH -J SF
#SBATCH --cpus-per-task=8
#SBATCH --output=/camcan/results/studyforrest/slurm-%j.out
#SBATCH --error=/camcan/results/studyforrest/slurm-%j.err

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python" 

echo "In flexiblehubs_1job parameters CHUNKLEN=$CHUNKLEN and SUBJ=$SUBJ"
${PYTHON} /camcan/FPNflexiblehubs/FPNflexiblehubs/flexiblehubsrunSF.py $CHUNKLEN $SUBJ
