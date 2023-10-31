# check which subsets are mismatched
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

expression_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_counts_expression_numpy_arrays"

sample_id_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_counts_sample_ID_arrays"

for num in range(0,718):
  sample_subset = pd.read_csv(sample_id_output_dir + "/refine.bio_subset_" + str(num)  + "_sampleIDs.txt", header = None)
  expression_subset = np.load(expression_output_dir + "/refine.bio_subset" + str(num) + "_counts_data.npy")
  if sample_subset.shape[0] != expression_subset.shape[0]:
    print("subset " + str(num) + " dimensions do not match")

print('The time it took this script to run is',time.time()-tic)

