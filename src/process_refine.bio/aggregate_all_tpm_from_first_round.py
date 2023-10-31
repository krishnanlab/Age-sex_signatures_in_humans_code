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

expression_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_tpm_expression_numpy_arrays"

expression_output = np.load(expression_output_dir + "/refine.bio_subset0_TPM_data.npy")

for num in range(1,718):
  tmp = np.load(expression_output_dir + "/refine.bio_subset" + str(num)  + "_TPM_data.npy")
  expression_output = np.concatenate((expression_output, tmp), axis = 0)
  print(str(num))


np.save(expression_output_dir + "/refine.bio_tpm_expression.npy", arr = expression_output)

print('The time it took this script to run is',time.time()-tic)



