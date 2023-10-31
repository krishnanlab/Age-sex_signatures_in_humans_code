import numpy as np
import pandas as pd 
from scipy.stats import zscore
from sklearn.metrics import balanced_accuracy_score
import matplotlib.pyplot as plt

##################READ DATA######################
#load gene expression data
gene_exp = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/X_2019-07-31.npy", allow_pickle=True)

#read in gene ids
gene_ids = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/X_GeneIDs_2019-07-31.npy", allow_pickle=True)

#read in sample ids
sample_ids = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/key_2019-07-31.npy", allow_pickle=True)
sample_ids = pd.DataFrame(sample_ids, columns=["gse", "gsm"])

#read in gene/chromosome data
gene_chr = pd.read_csv("/mnt/home/john3491/projects/age-sex-prediction/data/clean_mygeneinfo_homo_sapiens_entrez_genes_chromosome_locations.tsv", sep="\t", header=0)

#read in standard/sample labels
labels = pd.read_csv("/mnt/home/john3491/projects/age-sex-prediction/data/gpl570.gse-gsm-sex-age_group-age-tissue-disease.tsv", sep="\t", header=0)

############################GET Z SCORES/GENERATE STEPS##################
#calculate z scores in each column (gene)
gene_exp_zscores = np.apply_along_axis(zscore, 0, gene_exp, nan_policy='omit')

#find max and min
max_z = np.nanmax(gene_exp_zscores)
min_z = np.nanmin(gene_exp_zscores)

#generate steps from min to max z score by 0.1 to use for all genes
steps = np.arange(round(min_z, 1),round(max_z, 1), 0.2)

###################FILTER GENE EXP DATA, SAMPLES & GENES###################
#get X/Y genes from gene_chr, change to numpy array
gene_chr = gene_chr[(gene_chr['chromosome'] == 'X') | (gene_chr['chromosome'] == 'Y')]
genes = gene_chr['gene']
genes = genes.to_numpy()

#filter gene_exp_zscores to only include genes we want
#first makes gene_ids character array
gene_ids = gene_ids.astype('str')
#replace _at that is at the end of every gene for whatever reason
gene_ids = np.char.replace(gene_ids, '_at', '')
#get genes of interest which are in the gene_exp array:
#apparently it does not matter that gene_ids are str and genes are int64
gene_positions = np.isin(gene_ids, genes)
#subset gene_exp_zscores with gene_positions 
gene_exp_zscores = gene_exp_zscores[:, gene_positions]
#subset gene_ids to be same subset as gene_exp_zscores
gene_ids = gene_ids[gene_positions]

#subset samples since we have labels for like 19000 and may as well make it smaller
labeled_samples = labels['gsm']
sample_ids = sample_ids['gsm']
#get samples that have labels 
sample_positions = np.isin(sample_ids, labeled_samples)
#subset gene_exp_zscores to only have labeled samples
gene_exp_zscores = gene_exp_zscores[sample_positions, :]
#subset sample_ids to be same subset as gene_exp_zscores
sample_ids = sample_ids[sample_positions]

#get age groups
age_groups = labels.age_group.unique()
#get rid of samples that are in labels but not in gene_exp_zscores
labels = labels[labels['gsm'].isin(sample_ids)]
#create column for labels to give to balanced accuracy function
#female = T and male = F
labels['y_true'] = np.where(labels.sex == 'female', True, False)

final_gene_cuts = pd.DataFrame()

for ag in age_groups:
    #filter labels for given ag
    ag_labels = labels[labels['age_group']==ag]
    #get list of samples from ag
    ag_samples = ag_labels['gsm']
    #filter gene_exp_zscores for ag:
    #get boolean vec first
    ag_sample_positions = np.isin(sample_ids, ag_samples)
    #subset gene_exp_zscores to samples in ag
    ag_gene_expz = gene_exp_zscores[ag_sample_positions, :]
    #get sample_ids for ag_gene_expz
    ag_sample_ids = sample_ids[ag_sample_positions]

    #order labels with ag_sample_ids
    ag_labels = ag_labels.set_index('gsm')
    ag_labels = ag_labels.reindex(ag_sample_ids)
    #pull y_true for checking accuracy 
    ag_y_true = ag_labels['y_true']

    #get number cols in data
    col_num = ag_gene_expz.shape[1]

    #lists to append to
    best_cutoffs = []
    best_bas = []

    for col in np.arange(col_num):
        #get one whole column (gene)
        gene_values = ag_gene_expz[:, col]
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

        if gene_df.shape[0] == 1:
            best_cut = gene_df.iat[0,0]
        if gene_df.shape[0] > 1:
            best_cut = np.median(gene_df['cutoff'])
        
        best_cutoffs.append(best_cut)
        best_bas.append(high_ba)

    ag_gene_cuts = pd.DataFrame({'age_group':ag, 'gene':gene_ids, 'cutoff':best_cutoffs, 'balanced_accuracy':best_bas})
    final_gene_cuts = final_gene_cuts.append(ag_gene_cuts)

#done with entire ag loop here
#find best balanced_accuracies in each age group to plot
top15_gene_cuts = final_gene_cuts.groupby('age_group', group_keys=False).apply(pd.DataFrame.nlargest, n=15, columns='balanced_accuracy')

for ag in age_groups:
    #filter labels for given ag
    ag_labels = labels[labels['age_group']==ag]
    #get list of samples from ag
    ag_samples = ag_labels['gsm']
    #filter gene_exp_zscores for ag:
    #get boolean vec first
    ag_sample_positions = np.isin(sample_ids, ag_samples)
    #subset gene_exp_zscores to samples in ag
    ag_gene_expz = gene_exp_zscores[ag_sample_positions, :]
    #get sample_ids for ag_gene_expz
    ag_sample_ids = sample_ids[ag_sample_positions]

    #order labels with ag_sample_ids
    ag_labels = ag_labels.set_index('gsm')
    ag_labels = ag_labels.reindex(ag_sample_ids)
    #pull sex for plotting
    ag_sex = ag_labels['sex']

    #filter top15_gene_cuts to get only age group
    ag_top15_gene_cuts = top15_gene_cuts[top15_gene_cuts['age_group']==ag]
    #pull out top genes for plotting
    ag_best_genes = ag_top15_gene_cuts['gene']

    #make dictionary to append to
    ag_best_genes_exp_dic = {}

    for gene in ag_best_genes:
        #find gene in numpy array of gene exp
        #first find col number with gene_ids list
        gene_pos = np.where(gene_ids == gene)
        #pull gene exp column
        gene_values = ag_gene_expz[:, gene_pos]
        #get female data
        female_pos = np.where(ag_sex == 'female')
        female_exp_vals = gene_values[female_pos]
        #get male data
        male_pos = np.where(ag_sex == 'male')
        male_exp_vals = gene_values[male_pos]
        #append gene exp to dict
        ag_best_genes_exp_dic['female_' + gene] = female_exp_vals
        ag_best_genes_exp_dic['male_' + gene] = male_exp_vals
    
    #make figure
    fig, ax = plt.subplots(5,3, figsize = (8,11))
    fig.subplots_adjust(hspace=0.05, wspace=0.1)
    for i in range(5):
        for j in range(3):
            #get gene from ag_top15_gene_cuts
            gene = ag_top15_gene_cuts.iat[((i*3)+j), 1]
            #get balanced accuracy from ag_top15_gene_cuts (rounded to 2 dec places)
            ba = round(ag_top15_gene_cuts.iat[((i*3)+j), 3], 2)
            #get cutoff from ag_top15_gene_cuts (rounded to 2 dec places)
            cu = round(ag_top15_gene_cuts.iat[((i*3)+j), 2], 2)
            ax[i, j].hist([ag_best_genes_exp_dic['female_' + gene], ag_best_genes_exp_dic['male_' + gene]], alpha = 0.5, histtype = 'stepfilled', color = ['#ff7f0e', '#1f77b4'], weights = [np.ones(len(ag_best_genes_exp_dic['female_' + gene])) / len(ag_best_genes_exp_dic['female_' + gene]), np.ones(len(ag_best_genes_exp_dic['male_' + gene])) / len(ag_best_genes_exp_dic['male_' + gene])])
            ax[i, j].annotate(gene, xy=(0.05,0.9), xycoords = 'axes fraction')
            ax[i, j].annotate(ba, xy=(0.05,0.8), xycoords = 'axes fraction')
            ax[i, j].axvline(x = cu, color = 'black', linestyle = '--')
            ax[i, j].annotate(cu, xy=(cu, 0.4), xycoords = 'data')
        
    fig.suptitle("Top 15 genes for Naive Sex Prediction in the " + ag + " Age Group", y=1.01, fontsize=14)
    plt.tight_layout()
    plt.savefig('/mnt/home/john3491/projects/age-sex-prediction/results/naive_sex_prediction/' + ag + '_top_15_naive_sex_prediction_genes.pdf', format='pdf', dpi = 400, bbox_inches = 'tight')