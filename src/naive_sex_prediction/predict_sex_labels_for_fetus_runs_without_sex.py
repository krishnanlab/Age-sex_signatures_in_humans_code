import numpy as np
import pandas as pd 
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
import time

# script time
tic = time.time()

# genes with best cutoffs to separate male from female
best_genes = [8653, 9086, 6192, 7404, 8287, 83869, 7544, 7503, 8284, 9087, 22829, 90665, 6736, 246119, 159119]

# results output path
out_path = "../../results/naive_sex_prediction/"

# read in standard/sample labels
labels = pd.read_csv("../../data/labels/sample-filtered_manually_annotated_refine.bio_sample_labels.tsv", sep="\t", header=0, encoding = 'unicode_escape')
labels = labels[labels['age_group'].notnull()]
unsexed_fetus_labels = labels[labels['sex'].isnull()]
unsexed_fetus_runs = unsexed_fetus_labels['run']

# read in gene ids
gene_ids = np.load("../../data/expression/refine.bio_geneIDs.npy", allow_pickle=True)

#read in sample ids
sample_ids = pd.read_csv("../../data/expression/refine.bio_sample_IDs_tpm.txt", header = None)

# load gene expression data
gene_exp = np.load("../../data/expression/refine.bio_TPM_expression.npy", allow_pickle=True)
gene_exp = np.arcsinh(gene_exp)

# subset to genes being used for prediction
gene_positions = np.isin(gene_ids, best_genes)
#subset gene_exp with gene_positions 
gene_exp = gene_exp[:, gene_positions]
#subset gene_ids to be same subset as gene_exp
gene_ids = gene_ids[gene_positions]

# subset to fetus samples without sex
tpm_keep_mask = sample_ids[0].isin(unsexed_fetus_runs).tolist()
unlabeled_sample_ids = sample_ids[tpm_keep_mask]
unlabeled_exp = gene_exp[tpm_keep_mask, :]

# subset to labeled samples of all age groups
tpm_keep_mask = sample_ids[0].isin(labels['run']).tolist()
labeled_sample_ids = sample_ids[tpm_keep_mask]
labeled_exp = gene_exp[tpm_keep_mask, :]

# fit scaler to data of all age groups
std_scale = StandardScaler().fit(labeled_exp)

# subset to labeled samples of fetus group only
labels = labels[labels['fine_age_group'] == 'fetus']
tpm_keep_mask = sample_ids[0].isin(labels['run']).tolist()
labeled_sample_ids = sample_ids[tpm_keep_mask]
labeled_exp = gene_exp[tpm_keep_mask, :]

#get rid of samples that are in labels but not in gene_exp_zscores
labels = labels[labels['run'].isin(labeled_sample_ids[0])]

# z score fetus data using fit from all age groups
labeled_exp = std_scale.transform(labeled_exp)
unlabeled_exp = std_scale.transform(unlabeled_exp)

# create column for labels to give to balanced accuracy function
# female = T and male = F
labels['y_true'] = np.where(labels.sex == 'female', True, False)

# get sample_ids for labeled_exp
ag_sample_ids = labeled_sample_ids[0]
# order labels with ag_sample_ids
labels = labels.set_index('run')
labels = labels.reindex(ag_sample_ids)
# pull y_true for checking accuracy 
ag_y_true = labels['y_true']
# get number cols in data
col_num = labeled_exp.shape[1]

# find max and min
max_z = np.nanmax(labeled_exp)
min_z = np.nanmin(labeled_exp)
# generate steps from min to max z score by 0.1 to use for all genes
steps = np.arange(round(min_z, 1),round(max_z, 1), 0.1)

# lists to append to
best_cutoffs = []
best_bas = []
higher_sex = []

for col in np.arange(col_num):
    #get one whole column (gene)
    gene_values = labeled_exp[:, col]
    #find min/max values
    max_gene_z = np.nanmax(gene_values)
    min_gene_z = np.nanmin(gene_values)
    #subset step size vector to only cover values covered by gene
    gene_steps = steps[(steps > min_gene_z) & (steps < max_gene_z)]
    
    #make lists to append to
    gene_cutoffs = []
    gene_female_higher_ba = []
    gene_male_higher_ba = []
        
    if len(gene_steps) == 0:
        gene_cutoffs.append('none')
        gene_female_higher_ba.append(0)
        gene_male_higher_ba.append(0)
        
    if len(gene_steps) != 0:
        for cutoff in gene_steps:
            #get balanced accuracy if you classify everything above
            #cutoff as female
            pred = np.where(gene_values > cutoff, True, False)
            fh = balanced_accuracy_score(ag_y_true, pred)
            #get balanced accuracy if you classify everything above
            #cutoff as male
            pred = np.where(gene_values > cutoff, False, True)
            mh = balanced_accuracy_score(ag_y_true, pred)

            gene_cutoffs.append(cutoff)
            gene_female_higher_ba.append(fh)
            gene_male_higher_ba.append(mh)
        
    #put lists made in cutoff loop into one df
    gene_df = pd.DataFrame({'cutoff':gene_cutoffs, 'female_higher':gene_female_higher_ba, 'male_higher': gene_male_higher_ba})
    #find highest balanced accuracy
    high_ba = np.nanmax(gene_df[['female_higher','male_higher']])
    #get rows where highest_ba exists
    gene_df = gene_df[(gene_df['female_higher'] == high_ba) | (gene_df['male_higher'] == high_ba)]
        
    if gene_df['female_higher'].max() > gene_df['male_higher'].max():
        hs = 'female'
    if gene_df['male_higher'].max() > gene_df['female_higher'].max():
        hs = 'male'

    if gene_df.shape[0] == 1:
        best_cut = gene_df.iat[0,0]
    if gene_df.shape[0] > 1:
        best_cut = np.median(gene_df['cutoff'])
        
    best_cutoffs.append(best_cut)
    best_bas.append(high_ba)
    higher_sex.append(hs)

gene_cuts = pd.DataFrame({'age_group':'fetus', 'gene':gene_ids, 'cutoff':best_cutoffs, 'balanced_accuracy':best_bas, 'higher_sex':higher_sex})

gene_cuts.to_csv(out_path + "top_fetus_gene_cutoffs_and_balanced_accuracies.tsv", sep = "\t", header = True, index = False)

# output transformed expression
fetus_wo_sex_info = pd.DataFrame(unlabeled_exp, columns = gene_ids.tolist())
fetus_wo_sex_info['run'] = unlabeled_sample_ids[0].tolist()

fetus_wo_sex_info.to_csv(out_path + "std_scaled_fetus_samples_wo_sex_labels_expression_values.tsv", sep = "\t", header = True, index = False)

# script time
print('The time it took this script to run is',time.time()-tic)
