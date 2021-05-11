#!/bin/bash
#SBATCH -J FH-CC
#SBATCH --cpus-per-task=8
#SBATCH --output=/camcan/FPNflexiblehubs/FPNflexiblehubs/res/camcan/slurm-%j.out
#SBATCH --error=/camcan/FPNflexiblehubs/FPNflexiblehubs/res/camcan/slurm-%j.err

PYTHON="/camcan/anaconda3/envs/chiara_env/bin/python"


echo "In flexiblehubs_1jobPoolSubjsCC parameters CHUNKLEN=$CHUNKLEN"
${PYTHON} /camcan/FPNflexiblehubs/FPNflexiblehubs/flexiblehubsrunPoolSubjsCC.py $CHUNKLEN
