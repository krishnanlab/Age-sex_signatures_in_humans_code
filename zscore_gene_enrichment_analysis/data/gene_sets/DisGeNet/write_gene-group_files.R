library(tidyverse)

dfull <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_full.tsv")
dcur <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_curated.tsv")

dfull %>% 
  select(geneId, diseaseName) %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_full_gene-group_file.tsv")

dcur %>% 
  select(geneId, diseaseName) %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_curated_gene-group_file.tsv")
