#!/bin/bash -login
### define resources needed:
### Request the time
#SBATCH --time=28:00:00
### Request the total amount of memory
#SBATCH --mem=100GB
### Request the number of nodes
#SBATCH --nodes=1
### Request the number of cpus
#SBATCH --cpus-per-task=1
### Request to use our buyin nodes
#SBATCH --account=wang-krishnan

#this is for unzipping the file from refine.bio using wget
#change to current dir
cd $SLURM_SUBMIT_DIR

unzip HOMO_SAPIENS_2_1601762743.zip