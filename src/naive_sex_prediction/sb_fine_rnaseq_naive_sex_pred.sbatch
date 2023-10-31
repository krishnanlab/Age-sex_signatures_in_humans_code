#!/bin/sh -login
#SBATCH --time=04:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account=wang-krishnan

cd $SLURM_SUBMIT_DIR

source /mnt/home/john3491/.bashrc
export PATH="/mnt/home/john3491/anaconda3/bin:$PATH"

which python

python naive_sex_pred_by_age_group.py --data_type=rnaseq --group_type=fine
