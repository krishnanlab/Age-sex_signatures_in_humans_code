# this script aggregates the npy arrays from aggregate_tpm_data_first_round.py

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

sample_id_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_tpm_sample_ID_arrays"

sample_output = pd.read_csv(sample_id_output_dir + "/refine.bio_subset_0_sampleIDs.txt", header = None)

for num in range(1,718):
  tmp = pd.read_csv(sample_id_output_dir + "/refine.bio_subset_" + str(num)  + "_sampleIDs.txt", header = None)
  sample_output = pd.concat([sample_output, tmp], axis = 0)
  print(str(num))


sample_output.to_csv(sample_id_output_dir + "/refine.bio_tpm_sample_IDs.txt", header = False, index = False)

print('The time it took this script to run is',time.time()-tic)



