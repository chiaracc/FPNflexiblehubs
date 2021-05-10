#!/bin/bash
#SBATCH -J FH-SF
#SBATCH --cpus-per-task=8
#SBATCH --output=/camcan/FPNflexiblehubs/FPNflexiblehubs/res/studyforrest/slurm-%j.out
#SBATCH --error=/camcan/FPNflexiblehubs/FPNflexiblehubs/res/studyforrest/slurm-%j.err

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"


echo "In flexiblehubs_1jobPoolSubjsSF parameters CHUNKLEN=$CHUNKLEN"
${PYTHON} /camcan/FPNflexiblehubs/FPNflexiblehubs/flexiblehubsrunPoolSubjsSF.py $CHUNKLEN
