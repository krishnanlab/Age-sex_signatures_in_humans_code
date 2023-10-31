library(tidyverse)

gtex_tis_spec_scores <- read_tsv("/mnt/research/compbio/krishnanlab/data/GTEx/tissue-specific_genes/GTEx_expression_Entrez-tau_score-tissue-gene_zscore.tsv")
gtex_tis_spec_scores <- gtex_tis_spec_scores %>% 
  filter(tau_score > 0.8) %>% 
  filter(robust_zscore > 2) %>% 
  select(Entrez, tissue) %>% 
  write_tsv("/mnt/home/john3491/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/gtex_tissue-spec/gtex_tissue-spec_gene-group_file.tsv")


