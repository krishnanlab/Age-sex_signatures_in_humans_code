import numpy as np
import pandas as pd 
from scipy.stats import zscore
from sklearn.metrics import balanced_accuracy_score
import matplotlib.pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser()                                               
parser.add_argument("--data_type", help="rnaseq or microarrray", type=str, required=True)
args = parser.parse_args()

# script time
tic = time.time()

##################READ DATA######################
if args.data_type == 'microarray':
    #load gene expression data
    gene_exp = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/gpl570_gene_expression.npy", allow_pickle=True)
    #read in gene ids
    gene_ids = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/gene_IDs.npy", allow_pickle=True)

    #read in sample ids
    sample_ids = pd.read_csv("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/sample_IDs.txt", header = None)

    #read in standard/sample labels
    labels = pd.read_csv("/mnt/home/john3491/projects/age-sex-prediction/blood-only_analysis/data/blood_labels/gpl570_blood_samples.tsv", sep="\t", header=0)
    
if args.data_type == 'rnaseq':
    #load gene expression data
    gene_exp = np.load("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/data/refine.bio_tpm_expression_sample-filtered.npy", allow_pickle=True)
    #read in gene ids
    gene_ids = np.load("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/data/refine.bio_geneIDs.npy", allow_pickle=True)

    #read in sample ids
    sample_ids = pd.read_csv("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/data/refine.bio_filtered-sample_IDs_tpm.txt", header = None)
    
    #read in standard/sample labels
    labels = pd.read_csv("/mnt/home/john3491/projects/age-sex-prediction/blood-only_analysis/data/blood_labels/refine.bio_blood_samples.tsv", sep="\t", header=0, encoding = 'unicode_escape')
    labels = labels.rename(columns={"run": "gsm", "experiment": "gse"})
    labels = labels[labels['sex_inferred_wXY'] != True]

#subset samples
labeled_samples = labels['gsm']
#get samples that have labels 
sample_positions = np.isin(sample_ids[0], labeled_samples)
#subset gene_exp to only have labeled samples
gene_exp = gene_exp[sample_positions, :]
#get rid of second 1 dimension so next line works
sample_positions = sample_positions.squeeze()
#subset sample_ids to be same subset as gene_exp
sample_ids = sample_ids[sample_positions]

if args.data_type == 'rnaseq':
    gene_exp = np.arcsinh(gene_exp)
#calculate z scores in each column (gene)
gene_exp_zscores = np.apply_along_axis(zscore, 0, gene_exp, nan_policy='omit')
if args.data_type == 'rnaseq':
    #remove nan columns where not enough variation to do z scores
    no_nan_gene_positions = sum(np.isnan(gene_exp_zscores)) == 0
    gene_exp_zscores = gene_exp_zscores[:, no_nan_gene_positions]
    gene_ids = gene_ids[no_nan_gene_positions]

#get age groups
age_groups = ['adolescent', 'child', 'adult', 'young_adult', 'middle_adult', 'young_child', 'older_adult', 'old_adult', 'infant']
    
#create column for labels to give to balanced accuracy function
#female = T and male = F
labels['y_true'] = np.where(labels.sex == 'female', True, False)

final_gene_cuts = pd.DataFrame()

for ag in age_groups:
    #filter labels for given ag
    ag_labels = labels[labels['fine_age_group']==ag]
    #get list of samples from ag
    ag_samples = ag_labels['gsm']
    #filter gene_exp_zscores for ag:
    #get boolean vec first
    ag_sample_positions = np.isin(sample_ids[0], ag_samples)

    #subset gene_exp_zscores to samples in ag
    ag_gene_expz = gene_exp_zscores[ag_sample_positions, :]
    #get sample_ids for ag_gene_expz
    ag_sample_ids = sample_ids[ag_sample_positions]

    #order labels with ag_sample_ids
    ag_labels = ag_labels.set_index('gsm')
    ag_labels = ag_labels.reindex(ag_sample_ids[0])
    #pull y_true for checking accuracy 
    ag_y_true = ag_labels['y_true']

    #get number cols in data
    col_num = ag_gene_expz.shape[1]

    #lists to append to
    best_cutoffs = []
    best_bas = []
    higher_sex = []

    for col in np.arange(col_num):
        #get one whole column (gene)
        gene_values = ag_gene_expz[:, col]
        #find min/max values
        max_gene_z = np.nanmax(gene_values)
        min_gene_z = np.nanmin(gene_values)
        #subset step size vector to only cover values covered by gene
        gene_steps = np.linspace(start = round(min_gene_z, 2), stop = round(max_gene_z, 2), num = 20)
    
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

    ag_gene_cuts = pd.DataFrame({'age_group':ag, 'gene':gene_ids, 'cutoff':best_cutoffs, 'balanced_accuracy':best_bas, 'higher_sex':higher_sex})
    final_gene_cuts = final_gene_cuts.append(ag_gene_cuts)

final_gene_cuts.to_csv("/mnt/home/john3491/projects/age-sex-prediction/blood-only_analysis/data/naive_sex_prediction_signed_ranks/" + args.data_type + "_fine_age_groups_final_20-step_gene_cut_data_all_chromosomes.tsv", sep = "\t", header = True, index = False)

# script time
print('The time it took this script to run is',time.time()-tic)
