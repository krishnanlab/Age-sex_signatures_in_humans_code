#!/bin/bash -login
### define resources needed:
### Request the time 
#SBATCH --time=08:00:00
### Request the total amount of memory
#SBATCH --mem=50GB
### Request the number of nodes
#SBATCH --nodes=1
### Request the number of cpus 
#SBATCH --cpus-per-task=1
### Request to use our buyin nodes
#SBATCH --account=wang-krishnan

#change to current dir
cd $SLURM_SUBMIT_DIR

module purge
module load GCC/8.3.0 OpenMPI/3.1.4
module load R/4.0.2

Rscript find_jaccard_overlap_positive-negative_genes_age_models.R

