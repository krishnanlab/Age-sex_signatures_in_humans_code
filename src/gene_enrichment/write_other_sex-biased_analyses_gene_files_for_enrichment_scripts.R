
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

out_path <- "/mnt/home/john3491/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/"
#### read in data ####
# looking at X/Y genes I decided negative log fold changes
# are female-biased genes and positive are male-biased
sagd <- read_delim("/mnt/research/compbio/krishnanlab/data/SAGD/homo_sapiens/SAGD_human_sex-biased_entrez_genes_L2FCover2_padj-under0.05.tsv",
                   delim = "\t", col_names = T)

# positive effect size is female
gtex_genes <- read_delim("/mnt/research/compbio/krishnanlab/data/GTEx/sex-biased_analyses/GTEx_Analysis_v8_sbgenes/gtex_tissue_sex-biased_entrez_genes.tsv",
                         delim = "\t", col_names = T)
gtex_eqtls <- read_delim("/mnt/research/compbio/krishnanlab/data/GTEx/sex-biased_analyses/GTEx_Analysis_v8_sbeQTLs/gtex_tissue_sex-biased_entrez_eQTLs.tsv",
                         delim = "\t", col_names = T)

guo <- read_delim("/mnt/research/compbio/krishnanlab/data/sex-biased_Guo_et_al/guo_et_al_entrez_sex-baised_genes.tsv",
                  delim = "\t", col_names = T)

#### write sagd file ####
sagd$Tissue <- gsub(" ", "-", sagd$Tissue)
sagd <- sagd %>% 
  mutate(bias = ifelse(log2FoldChange < 0, "female", "male")) %>% 
  mutate(group = paste(Tissue, Stage, bias, sep = "_")) %>% 
  select(entrez, group) %>% 
  distinct()

sagd %>% write_delim(paste0(out_path, "sagd_gene-group_file.tsv"),
                     delim = "\t", col_names = T)

#### write gtex files ####
gtex_genes$tissue <- gsub("_", "-", gtex_genes$tissue)
gtex_genes <- gtex_genes %>% 
  mutate(bias = ifelse(effsize > 0, "female", "male")) %>% 
  mutate(group = paste(tissue, bias, sep = "_")) %>% 
  select(entrez, group) %>% 
  distinct()

gtex_genes %>% write_delim(paste0(out_path, "gtex-genes_gene-group_file.tsv"),
                           delim = "\t", col_names = T)

#### write guo files ####
guo$tissue <- gsub(" ", "-", guo$tissue)
guo$sex_bias <- gsub("F", "female", guo$sex_bias)
guo$sex_bias <- gsub("M", "male", guo$sex_bias)

guo <- guo %>% 
  mutate(group = paste(tissue, sex_bias, sep = "_")) %>% 
  select(entrez, group) %>% 
  distinct()

guo %>% write_delim(paste0(out_path, "guo_gene-group_file.tsv"),
                    delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
