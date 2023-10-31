# this script writes samples that are not over half zero in counts
# these are samples I will keep for the age/sex project

import pickle
import numpy as np
import datetime
import sys
import time
import pandas as pd
import os
import re

# script time
tic = time.time()

path = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/data/"

######### find samples that pass filter #############
# samples are rows, genes are columns
counts = np.load(path + "refine.bio_counts_expression.npy")

# read in sample IDs for counts
counts_samples = pd.read_csv(path + "refine.bio_counts_sample_IDs.txt", sep = "\t", header = None)

# find samples > 50% zero
# array of the sum of nonzero elements in each row
row_nonzeros = np.count_nonzero(counts, axis = 1)
# find indices of samples over half zero (25570 genes)
keep_mask = row_nonzeros > 12784
# subset counts samples to samples that pass filter
pass_samples = counts_samples[keep_mask]

# write samples to keep
pass_samples.to_csv("/mnt/home/john3491/projects/age-sex-prediction/data/refine.bio/refine.bio_samples_over_half_nonzero_counts.txt", header = False, index = False)

############ read in tpm data/samples to subset ##############
tpm = np.load(path + "refine.bio_tpm_expression.npy")

# sample IDs for tpm
tpm_samples = pd.read_csv(path + "refine.bio_tpm_sample_IDs.txt", sep = "\t", header = None)

# subset tpm samples and expression to pass samples
tpm_keep_mask = tpm_samples[0].isin(pass_samples[0]).tolist()
tpm_samples = tpm_samples[tpm_keep_mask]
tpm = tpm[tpm_keep_mask, :]

############### filter out single cell and cell line projects ###############
home_path = "/mnt/home/john3491/projects/age-sex-prediction/data/"
bad_projects = pd.read_csv(home_path + "labels/bad_refine.bio_projects_to_filter_out.txt", sep = "\t", header = None)
bad_projects = bad_projects[0].tolist()
refinebio_metadata = pd.read_csv(home_path + "refine.bio/full_refine.bio-metasra_aggregated_metadata.tsv", sep ="\t", low_memory=False)

bad_runs = refinebio_metadata[refinebio_metadata.experiment_accession.isin(bad_projects)]
bad_runs = bad_runs['run'].tolist()

keep_mask = tpm_samples[0].isin(bad_runs) == False
keep_mask = keep_mask.tolist()

# subset tpm samples and expression to pass samples
tpm_samples = tpm_samples[keep_mask]
tpm = tpm[keep_mask, :]

# write sample filtered tpm and samples
tpm_samples.to_csv(path + "refine.bio_filtered-sample_IDs_tpm.txt", header = False, index = False)
np.save(path + "refine.bio_tpm_expression_sample-filtered.npy", arr = tpm)

# script time
print('The time it took this script to run is',time.time()-tic)
