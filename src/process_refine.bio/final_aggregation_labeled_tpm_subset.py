# put together the subsections of tpm data
# from aggregate_labeled_subset_tpm_data.py

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

# refine.bio tpm dir
data_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/tpm"

# sample list
samples = []

for i in list(range(0,28)):
  with open(data_dir + "/labeled_subset_" + str(i) + "_sampleIDs.txt") as f:
    lines = f.read().splitlines()
  samples.extend(lines)
  
tpm = np.load(data_dir + "/refine.bio_labeled_subset0_TPM_data.npy")  
for i in list(range(1,28)):
  subset = np.load(data_dir + "/refine.bio_labeled_subset" + str(i) + "_TPM_data.npy")
  tpm = np.concatenate((tpm, subset), axis=0)
  
# write to file
samples = pd.DataFrame({'samples': samples})
samples.to_csv(data_dir + "/labeled_subset_sampleIDs.txt", header = False, index = False)

np.save(data_dir + "/refine.bio_labeled_subset_TPM_data.npy", arr = tpm)

# print shapes
print("the length of sample list is", len(samples))
print("the shape of the tpm data is", tpm.shape)

# script time
print('The time it took this script to run is',time.time()-tic)
