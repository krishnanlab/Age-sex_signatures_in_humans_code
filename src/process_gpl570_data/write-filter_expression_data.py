import numpy as np
import pandas as pd 

#load gene expression data
gene_exp = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/X_2019-07-31.npy", allow_pickle=True)
tmp = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/y_2019-07-31.npy", allow_pickle=True)

gene_exp = np.concatenate((gene_exp, tmp), axis = 1)

#read in sample ids
sample_ids = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/key_2019-07-31.npy", allow_pickle=True)
sample_ids = pd.DataFrame(sample_ids, columns=["gse", "gsm"])

# filter out nan samples (GSM1299523, GSM1274714, )
nan_positions = np.isnan(gene_exp[:, 0])
gene_exp = gene_exp[nan_positions == False, :]
sample_ids = sample_ids[~nan_positions]

# save filtered data
np.save("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/sample_info.npy", arr=sample_ids)
np.save("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/gpl570_gene_expression.npy", arr=gene_exp)

sample_ids = sample_ids['gsm']
np.save("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/sample_IDs.npy", arr=sample_ids)
sample_ids.to_csv("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/sample_IDs.txt", sep = "\t", header = None, index = None)

#read, combine, save gene IDs too
gene_ids = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/X_GeneIDs_2019-07-31.npy", allow_pickle=True)
tmp = np.load("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/y_GeneIDs_2019-07-31.npy", allow_pickle=True)

gene_ids = np.concatenate((gene_ids, tmp))
gene_ids = gene_ids.tolist()
gene_ids = [g.replace('_at', '') for g in gene_ids]

np.save("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/gene_IDs.npy", arr = gene_ids)

gene_ids = pd.DataFrame(gene_ids)
gene_ids.to_csv("/mnt/research/compbio/krishnanlab/data/GEO/2019-07-29_downloaded-files/age-sex_project/gene_IDs.txt", sep="\t", header=None, index=None)




