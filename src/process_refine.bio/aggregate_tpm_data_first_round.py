# this script aggregates the unzipped refine.bio data
# where each directory is a project and each file
# in a directory is a sample or the experiment metadata
# Need to convert ENST to Entrez genes as well

import pickle
import numpy as np
import datetime
import sys
import time
import pandas as pd
import os
import re
import argparse
from pandas.errors import ParserError

parser = argparse.ArgumentParser()                                               
parser.add_argument("section", help="section of samples to copy", type=int)
args = parser.parse_args()

# script time
tic = time.time()

# Entrez/ENST conversion
# which I am going to use to sum both
# tpm and counts data to gene level
# rename is to better merge with refine.bio data
new_file_info = []
with open('/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__ENST-to-Entrez__All-Mappings.tsv','r') as f:
    for idx, line in enumerate(f):
        col1id = line.strip().split('\t')[0]
        col2_list = line.strip().split('\t')[1].split(', ')
        for item in col2_list:
            new_file_info.append([col1id,item])

tpm = pd.DataFrame(new_file_info,columns=['Name','gene'])

# refine.bio projects dir
data_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/"

expression_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_tpm_expression_numpy_arrays"

sample_id_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/aggregated_tpm_sample_ID_arrays"
gene_id_output_dir = "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/tpm_gene_IDs"

if not os.path.isdir(expression_output_dir):
  os.makedirs(expression_output_dir)
  
if not os.path.isdir(sample_id_output_dir):
  os.makedirs(sample_id_output_dir)
  
if not os.path.isdir(gene_id_output_dir):
  os.makedirs(gene_id_output_dir)
  
# read list of dirs
experiments = pd.read_csv("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/all_refine.bio_experiment_directories.txt", sep = "\t", header = None)

# subset to correct section for job
ss = 9 * args.section
se = ss + 9
experiments = experiments[0][ss:se]

# sample list
samples = []

for exp in experiments:
  dir = data_dir + exp
  print(dir)
  files = os.listdir(dir)
  files = [i for i in files if i.endswith('.sf')]
  for file in files:
    # read in data: Name (of transcript) | Length | EffectiveLength | TPM | NumReads
    try: 
      dat = pd.read_csv(dir + "/" + file, sep = "\t")
    except ParserError:
      print(file + " could not be parsed")
      continue
    if list(dat.columns) != ['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads']:
      print(file + " does not have expected keys")
      continue
    # get sample ID from file name
    f = file.replace('_quant.sf', '')
    # merge tpm datatpm
    tdat = dat[['Name', 'TPM']]
    if tdat.TPM.dtype == 'O':
      print(file + " has nuisance numbers")
      continue
    tdat = tdat.rename(columns = {'TPM': f})
    tpm = tpm.merge(tdat, on = 'Name', how = 'left')
    # add to samples list
    samples.append(f)

print("finished loop in ",time.time()-tic)
# drop transcript column
tpm.drop('Name', axis = 1, inplace = True)
if list(tpm.columns.values)[1:] != samples:
  print("samples added incorrectly")
  
# summarize to gene level
tpm = tpm.groupby('gene').sum()
# reset the index because grouping makes 'gene' the index
tpm.reset_index(inplace = True)

print("finished summing to gene level in ", time.time()-tic) 

# write genes
genes = tpm['gene']
genes.to_csv(gene_id_output_dir + "/refine.bio_subset_" + str(args.section) + "_geneIDs.txt", header = False, index = False)

#write samples
samples = pd.DataFrame({'samples': samples})
samples.to_csv(sample_id_output_dir + "/refine.bio_subset_" + str(args.section) + "_sampleIDs.txt", header = False, index = False)

# drop gene column
tpm.drop('gene', axis = 1, inplace=True)

# convert to numpy array
tpm = tpm.to_numpy()

# transpose so samples are rows and genes are columns
tpm = tpm.transpose()

print("finished transpose in ",time.time()-tic) 

# save as numpy array
np.save(expression_output_dir + "/refine.bio_subset" + str(args.section) + "_TPM_data.npy", arr = tpm)

print('The time it took this script to run is',time.time()-tic)

