library(tidyverse)

gobps <- read_tsv("/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/Gene_Ontology/propagated_annotations/GO__propagated-annotations__Homo_sapiens__Entrez__BP__EXP_IDA_IPI_IMP_IGI_TAS_IC.tsv")

gobps <- gobps %>% 
  filter(gene_count > 9) %>% 
  filter(gene_count < 201)

gobps %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GOBPs/GO__propagated-annotations__term-size10-200__Homo_sapiens__Entrez__BP__EXP_IDA_IPI_IMP_IGI_TAS_IC.tsv",
            col_names = T)

gobps %>% 
  select(gene_id, go_id) %>% 
  separate_rows(gene_id, sep = ", ", convert = T) %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GOBPs/GOBP_gene-group_file.tsv",
            col_names = T)

