import numpy as np
import pandas as pd 
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.linear_model import LogisticRegression
import time
import argparse
import os

parser = argparse.ArgumentParser()                                               
parser.add_argument("--sex", help="male or female", type=str, required=True)
parser.add_argument("--pos_age_group", help="age group to train model for, see age_groups variable in script for options", type=str, required=True)
parser.add_argument("--data_type", help="rnaseq or microaarray", type=str, required=True)
parser.add_argument("--std_scale", help="with or without", type=str, required=True)
parser.add_argument('--asinh', dest='feature', action='store_true')
parser.add_argument('--no-asinh', dest='feature', action='store_false')
parser.set_defaults(feature=True)
args = parser.parse_args()

# script time
tic = time.time()

print("sex: " + args.sex)
print("pos_age_group: " + args.pos_age_group)
print("data_type: " + args.data_type)

# age_groups
age_groups = ["fetus","infant","young_child","child","adolescent","young_adult","adult","middle_adult","older_adult","old_adult","elderly"]
  
# data type
if args.data_type == "rnaseq":
  label_path = "../../data/folds/refine.bio_3F_CV_folds.tsv"
  expression_path = "../../data/expression/refine.bio_TPM_expression.npy"
  sample_path = "../../data/expression/refine.bio_sample_IDs_tpm.txt"
  gene_path = "../../data/expression/refine.bio_geneIDs.npy"

if args.data_type == "microarray":
  label_path = "../../data/folds/gpl570_3F_CV_folds.tsv"
  expression_path = "../../data/expression/gpl570_expression.npy"
  sample_path = "../../data/expression/gpl570_sample_IDs.txt"
  gene_path = "../../data/expression/gpl570_gene_IDs.npy"
  
# read in labels and folds
labels = pd.read_csv(label_path, sep = "\t")
if args.data_type == "microarray":
  labels.rename(columns={'gsm': 'run'}, inplace=True)

# load gene expression data
gene_exp = np.load(expression_path, allow_pickle=True)

if args.data_type == "rnaseq":
  if args.feature:
    gene_exp = np.arcsinh(gene_exp)

# read in sample ids
sample_ids = pd.read_csv(sample_path, header = None, sep = "\t")

# read in gene ids
gene_ids = np.load(gene_path)

# read in common genes between microarray and rnaseq
common_genes = pd.read_csv("../../data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", header = None)

# subset to common genes
common_gene_positions = np.isin(gene_ids, common_genes[0])
gene_ids = gene_ids[common_gene_positions]
gene_exp = gene_exp[:, common_gene_positions]


# subset to remove samples that don't pass filter
filter_mask = np.isin(labels['run'], sample_ids[0])
labels = labels[filter_mask]

# subset to correct sex
sex_labels = labels[labels['sex'] == args.sex]

nfolds = [1, 2, 3]

auroc_results = []
model_weights = []
positives = []
sample_probs_df_list = []

for fold in nfolds:
  tst_labels = sex_labels[sex_labels['fold'] == fold]
  trn_labels = sex_labels[sex_labels['fold'] != fold]
  
  tst_samples_mask = np.isin(sample_ids[0], tst_labels['run'])
  tst_sample_ids = sample_ids[tst_samples_mask]
  tst_data = gene_exp[tst_samples_mask, :]
  
  trn_samples_mask = np.isin(sample_ids[0], trn_labels['run'])
  trn_sample_ids = sample_ids[trn_samples_mask]
  trn_data = gene_exp[trn_samples_mask, :]
  
  # reindex labels to match sample order
  tst_labels = tst_labels.set_index('run')
  tst_labels = tst_labels.reindex(tst_sample_ids[0].tolist())
  
  trn_labels = trn_labels.set_index('run')
  trn_labels = trn_labels.reindex(trn_sample_ids[0].tolist())
  
  # add y_true column to labels
  tst_labels['y_true'] = np.where(tst_labels.age_group == args.pos_age_group, True, False).tolist()
  trn_labels['y_true'] = np.where(trn_labels.age_group == args.pos_age_group, True, False).tolist()
  
  # get vector of y labels for each
  y_tst = tst_labels.y_true.values
  y_trn = trn_labels.y_true.values
  
  # get number of positives
  npos = sum(y_tst)
  positives.append(npos)
  
  # scale train data
  if args.std_scale == "with":
    std_scale = StandardScaler().fit(trn_data)
    trn_data = std_scale.transform(trn_data)
    tst_data = std_scale.transform(tst_data)
  
  # train a model
  clf = LogisticRegression(penalty = 'elasticnet', max_iter = 30000, solver = 'saga', l1_ratio = 0.5)
  clf.fit(trn_data, y_trn)
  
  # save model weights
  weights = clf.coef_.flatten()
  # probability sample is a positive
  probs = clf.predict_proba(tst_data)[:,1]
  # auroc results
  alt_ag_aurocs = []
  for ag in age_groups:
    tst_lab = np.where(tst_labels.age_group == ag, True, False).tolist()
    auroc = roc_auc_score(tst_lab, probs)
    alt_ag_aurocs.append(auroc)
  
  model_weights.append(weights)
  auroc_results.append(alt_ag_aurocs)
  # sample probs
  tmp_sp = pd.DataFrame({'n_positives': npos, 'tst_fold': fold, 'pos_age_group': args.pos_age_group, 'sex': args.sex, 'y_true': y_tst, 'sample_probability': probs, 'run': tst_labels.index})
  sample_probs_df_list.append(tmp_sp)

# write results  
results_path = "../../results/age_prediction_by_sex/"+ args.data_type + "/output_files/"
if not os.path.exists(results_path):
   os.makedirs(results_path)

model_weights = pd.DataFrame(model_weights, columns = gene_ids.tolist())
model_weights['n_positives'] = positives
auroc_results = pd.DataFrame(auroc_results, columns = age_groups)
auroc_results['n_positives'] = positives
sample_probs = pd.concat(sample_probs_df_list, ignore_index=True)

if args.data_type == "microarray":
  auroc_results.to_csv(results_path + "3FCV_elasticnet_LR_auroc_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)
  model_weights.to_csv(results_path + "3FCV_elasticnet_LR_model_weights_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_weights", sep = "\t", header = True, index = False)
  sample_probs.to_csv(results_path + "3FCV_elasticnet_LR_sample_probabilities_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)

if args.data_type == "rnaseq":
  if args.feature: 
    auroc_results.to_csv(results_path + "3FCV_elasticnet_LR_auroc_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)
    model_weights.to_csv(results_path + "3FCV_elasticnet_LR_model_weights_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_weights.tsv", sep = "\t", header = True, index = False)
    sample_probs.to_csv(results_path + "3FCV_elasticnet_LR_sample_probabilities_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)
  if not args.feature:
    auroc_results.to_csv(results_path + "3FCV_elasticnet_LR_auroc_no-asinh_transformation_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)
    model_weights.to_csv(results_path + "3FCV_elasticnet_LR_model_weights_no-asinh_transformation_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_weights.tsv", sep = "\t", header = True, index = False)
    sample_probs.to_csv(results_path + "3FCV_elasticnet_LR_sample_probabilities_no-asinh_transformation_" + args.std_scale + "_std_scaling_" + args.sex + "_" + args.pos_age_group + "_results.tsv", sep = "\t", header = True, index = False)

# script time
print('The time it took this script to run is',time.time()-tic)
