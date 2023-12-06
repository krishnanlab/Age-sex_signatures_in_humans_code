import numpy as np
import pandas as pd 
from sklearn.metrics import average_precision_score
import time
# script time
tic = time.time()

fine_age_groups = ["fetus", "infant", "young_child", "child", "adolescent", "young_adult", "adult", "middle_adult", "older_adult", "old_adult", "elderly"]

rflabel_path = "../../data/folds/refine.bio_3F_CV_folds.tsv"
mflabel_path = "../../data/folds/gpl570_3F_CV_folds.tsv"

# read in labels
rflabels = pd.read_csv(rflabel_path, sep = "\t")
mflabels = pd.read_csv(mflabel_path, sep = "\t")

# adjust cols
rflabels.rename(columns = {'age_group':'run_fine_age_group'}, inplace = True)
mflabels.rename(columns = {'age_group':'run_fine_age_group', 'gsm':'run'}, inplace = True)

rflabels = rflabels[['run', 'run_fine_age_group']]
mflabels = mflabels[['run', 'run_fine_age_group']]

rsample_probs_path = "../../results/age_prediction_by_sex/rnaseq/sample_probabilites_from_model_with_chosen_parameters.tsv"
msample_probs_path = "../../results/age_prediction_by_sex/microarray/sample_probabilites_from_model_with_chosen_parameters.tsv"

rout_path = "../../results/age_prediction_by_sex/rnaseq/"
mout_path = "../../results/age_prediction_by_sex/microarray/"

# read in sample probs
rsample_probs = pd.read_csv(rsample_probs_path, sep = "\t")
msample_probs = pd.read_csv(msample_probs_path, sep = "\t")

# drop this col bc exactly the same as 'pos_age_group' col
rsample_probs.drop(columns=['age_group'], inplace = True)
msample_probs.drop(columns=['age_group'], inplace = True)

rfsample_probs = rsample_probs.merge(rflabels, how = 'left', on = 'run')

mfsample_probs = msample_probs.merge(mflabels, how = 'left', on = 'run')

rfsample_probs = rfsample_probs[rfsample_probs['group_type'] == 'fine']

mfsample_probs = mfsample_probs[mfsample_probs['group_type'] == 'fine']

##### rnaseq fine auprc results #####
auprc_results = pd.DataFrame({"fetus":[0.22, 0.22], "infant": [0.22, 0.22], "young_child": [0.22, 0.22], "child": [0.22, 0.22], "adolescent": [0.22, 0.22], "young_adult": [0.22, 0.22], "adult": [0.22, 0.22], "middle_adult": [0.22, 0.22], "older_adult": [0.22, 0.22], "old_adult": [0.22, 0.22], "elderly": [0.22, 0.22]})
all_models = rfsample_probs.parameters.unique()
for model in all_models:
  model_sample_probs = rfsample_probs[rfsample_probs['parameters'] == model]
  for fold in [1, 2, 3]:
    fold_df = model_sample_probs[model_sample_probs['tst_fold'] == fold]
    alt_ag_auprcs = []
    for ag in fine_age_groups:
      tst_lab = np.where(fold_df.run_fine_age_group == ag, 1, 0)
      auprc = average_precision_score(tst_lab, fold_df.sample_probability.values)
      alt_ag_auprcs.append(auprc)
    alt_ag_auprcs_df = pd.DataFrame([alt_ag_auprcs], columns = ["fetus", "infant", "young_child", "child", "adolescent", "young_adult", "adult", "middle_adult", "older_adult", "old_adult", "elderly"])
    auprc_results = pd.concat([auprc_results, alt_ag_auprcs_df])
 
auprc_results = auprc_results.iloc[2: , :]  
auprc_results['parameters'] = np.repeat(all_models, 3).tolist()
auprc_results['tst_fold'] = [1, 2, 3] * 22

auprc_results.to_csv(rout_path + "auprc_results.tsv", sep = "\t", header = True, index = False)

##### microarray fine auprc results #####
auprc_results = pd.DataFrame({"fetus":[0.22, 0.22], "infant": [0.22, 0.22], "young_child": [0.22, 0.22], "child": [0.22, 0.22], "adolescent": [0.22, 0.22], "young_adult": [0.22, 0.22], "adult": [0.22, 0.22], "middle_adult": [0.22, 0.22], "older_adult": [0.22, 0.22], "old_adult": [0.22, 0.22], "elderly": [0.22, 0.22]})
all_models = mfsample_probs.parameters.unique()
for model in all_models:
  model_sample_probs = mfsample_probs[mfsample_probs['parameters'] == model]
  for fold in [1, 2, 3]:
    fold_df = model_sample_probs[model_sample_probs['tst_fold'] == fold]
    alt_ag_auprcs = []
    for ag in fine_age_groups:
      tst_lab = np.where(fold_df.run_fine_age_group == ag, 1, 0)
      auprc = average_precision_score(tst_lab, fold_df.sample_probability.values)
      alt_ag_auprcs.append(auprc)
    alt_ag_auprcs_df = pd.DataFrame([alt_ag_auprcs], columns = ["fetus", "infant", "young_child", "child", "adolescent", "young_adult", "adult", "middle_adult", "older_adult", "old_adult", "elderly"])
    auprc_results = pd.concat([auprc_results, alt_ag_auprcs_df])
 
auprc_results = auprc_results.iloc[2: , :]  
auprc_results['parameters'] = np.repeat(all_models, 3).tolist()
auprc_results['tst_fold'] = [1, 2, 3] * 22

auprc_results.to_csv(mout_path + "auprc_results.tsv", sep = "\t", header = True, index = False)

# script time
print('The time it took this script to run is',time.time()-tic)
