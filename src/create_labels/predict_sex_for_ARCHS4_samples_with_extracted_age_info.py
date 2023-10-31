# this script gives a sex label to all samples in 
# human_ARCHS4_text_parsed_gse-gsm-sex-age-age_group.tsv that had
#age information extracted but not sex
import numpy as np
import pandas as pd 
from sklearn.metrics import balanced_accuracy_score

path = "/mnt/research/compbio/krishnanlab/data/rnaseq/archs4/age_sex_project/filtered_samples_over_half_zero_exp/transformed_data/"
#read in labels, get subset with sex
labels = pd.read_csv(path + "human_ARCHS4_text_parsed_gse-gsm-sex-age-age_group.tsv", sep="\t", header=0)
sex_labels = labels[labels['sex'].notnull()]
#read in sample info, get GSMs of gene_exp
sample_ids = pd.read_csv(path + "Sample_info.tsv", sep="\t", header=0)
sample_ids = sample_ids['GSM']
#read in gene ids of gene_exp
gene_ids = pd.read_csv(path + "human_GeneIDs.txt", sep="\t", header=None)
gene_ids = gene_ids[0]

gene_exp = np.load(path + "asinh_z-scored_expression_data.npy", allow_pickle=True)

#read in gene/chromosome data
gene_chr = pd.read_csv("/mnt/home/john3491/projects/age-sex-prediction/data/clean_mygeneinfo_homo_sapiens_entrez_genes_chromosome_locations.tsv", sep="\t", header=0)


#get GSMs with sex labels and use them to subset
sexed_gsms = sex_labels['gsm']
sexed_sample_positions = np.isin(sample_ids, sexed_gsms)
#subset gene_exp to only have sexed samples
sexed_gene_exp = gene_exp[sexed_sample_positions, :]
#subset sample_ids to be same subset as sexed_gene_exp
sexed_sample_ids = sample_ids[sexed_sample_positions]

#subset sexed gene exp to X/Y genes
#get X/Y genes from gene_chr
gene_chr = gene_chr[(gene_chr['chromosome'] == 'X') | (gene_chr['chromosome'] == 'Y')]
xygenes = gene_chr['gene']
gene_positions = np.isin(gene_ids, xygenes)
#subset sexed_gene_exp with gene_positions 
sexed_gene_exp = sexed_gene_exp[:, gene_positions]
#subset gene_ids to be same subset as sexed_gene_exp
xy_gene_ids = gene_ids[gene_positions]


#find max and min
max_z = np.nanmax(gene_exp)
min_z = np.nanmin(gene_exp)
#generate steps from min to max z score by 0.1 to use for all genes
steps = np.arange(round(min_z, 1),round(max_z, 1), 0.2)
#get number of cols (genes) to step through
col_num = sexed_gene_exp.shape[1]

#reorder sex labels to match sample IDs
sex_labels = sex_labels.set_index('gsm')
sex_labels = sex_labels.reindex(sexed_sample_ids)
#create column for labels to give to balanced accuracy function
#female = T and male = F
sex_labels['y_true'] = np.where(sex_labels.sex == 'female', True, False)

#lists to append to
best_cutoffs = []
best_bas = []
higher_sex = []

for col in np.arange(col_num):
    #get one whole column (gene)
    gene_values = sexed_gene_exp[:, col]
    #find min/max values
    max_gene_z = np.nanmax(gene_values)
    min_gene_z = np.nanmin(gene_values)
    #subset step size vector to only cover values covered by gene
    gene_steps = steps[(steps > min_gene_z) & (steps < max_gene_z)]
    #make lists to append to
    gene_cutoffs = []
    gene_female_higher_ba = []
    gene_male_higher_ba = []
    
    for cutoff in gene_steps:
        #get balanced accuracy if you classify everything above
        #cutoff as female
        pred = np.where(gene_values > cutoff, True, False)
        fh = balanced_accuracy_score(sex_labels['y_true'], pred)
        #get balanced accuracy if you classify everything above
        #cutoff as male
        pred = np.where(gene_values > cutoff, False, True)
        mh = balanced_accuracy_score(sex_labels['y_true'], pred)

        gene_cutoffs.append(cutoff)
        gene_female_higher_ba.append(fh)
        gene_male_higher_ba.append(mh)

    #put lists made in cutoff loop into one df
    gene_df = pd.DataFrame({'cutoff':gene_cutoffs, 'female_higher':gene_female_higher_ba, 'male_higher': gene_male_higher_ba})
    #find highest balanced accuracy
    high_ba = np.nanmax(gene_df[['female_higher','male_higher']])
    #get rows where highest_ba exists
    gene_df = gene_df[(gene_df['female_higher'] == high_ba) | (gene_df['male_higher'] == high_ba)]
    #pull out best cut
    if gene_df.shape[0] == 1:
        best_cut = gene_df.iat[0,0]
    if gene_df.shape[0] > 1:
        best_cut = np.median(gene_df['cutoff'])
    #pull out whether male or female expresion is higher at best cut
    if gene_df.iat[0,1] > gene_df.iat[0,2]:
        hs = 'female'
    if gene_df.iat[0,1] < gene_df.iat[0,2]:
        hs = 'male'
    
    best_cutoffs.append(best_cut)
    best_bas.append(high_ba)
    higher_sex.append(hs)
    
sex_cuts = pd.DataFrame({'gene': xy_gene_ids, 'cutoff':best_cutoffs, 'balanced_accuracy':best_bas, 'higher_sex':higher_sex})
    
high_sex_cuts = sex_cuts[sex_cuts['balanced_accuracy'] > 0.85]
    
##########Use genes found above to sex the rest of the samples with age group labels###########
unsexed_labels = labels[labels['sex'].isnull()]

#get GSMs with sex labels and use them to subset
unsexed_gsms = unsexed_labels['gsm']
unsexed_sample_positions = np.isin(sample_ids, unsexed_gsms)
#subset gene_exp to only have sexed samples
unsexed_gene_exp = gene_exp[unsexed_sample_positions, :]
#subset sample_ids to be same subset as sexed_gene_exp
unsexed_sample_ids = sample_ids[unsexed_sample_positions]

#subset data to genes with high balanced accuracy
high_gene_positions = np.isin(gene_ids, high_sex_cuts['gene'])
#subset sexed_gene_exp with gene_positions 
unsexed_gene_exp = unsexed_gene_exp[:, high_gene_positions]
#subset gene_ids to be same subset as sexed_gene_exp
high_gene_ids = gene_ids[high_gene_positions]

unsexed_labels = unsexed_labels.set_index('gsm')
unsexed_labels = unsexed_labels.reindex(unsexed_sample_ids)

labels_by_gene = unsexed_gene_exp

# IMPORTANT NOTE PLEASE READ
# This loop only works because I know all of the best genes cuts have
# higher expression for males
high_gene_num = unsexed_gene_exp.shape[1]
for i in np.arange(high_gene_num):
    #female = 1, male = 0
    labels_by_gene[:,i] = np.where(labels_by_gene[:,i] > high_sex_cuts.iat[i, 1], 0, 1)
    
labels_by_gene = labels_by_gene.sum(axis=1)
#sex column will now be 0-7, depending on how many genes declared it
# male or female
# sum of zero means it was 'male' by all 7 gene cutoffs while
# sum of 7 means it was 'female' by all 7 gene cutoffs
unsexed_labels['sex'] = labels_by_gene

unsexed_labels = unsexed_labels.reset_index()
unsexed_labels.to_csv(path + "samples_not_sexed_by_metadata_parsing.tsv", sep = "\t", header=True, index=False)


    

    
