library(readxl)
library(tidyverse)

sc <- read_excel("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sctype/ScTypeDB_full.xlsx") %>% 
  select(-geneSymbolmore2, -shortName) %>% 
  separate_rows(geneSymbolmore1, sep = ",")

gene_sym <- read_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/20step_naive_sex_prediction/rnaseq_signed_balanced_accuracy_data.tsv",
                     col_names = T) %>% 
  select(gene, symbol) %>% 
  distinct()

# choose one symbol for duplicate genes with diff symbols
gene_sym <- gene_sym %>% 
  filter(!(gene == 1394 & symbol == "LINC02210-CRHR1")) %>% 
  filter(!(gene == 23596 & symbol == "CHML")) %>% 
  filter(!(gene == 3117 & symbol == "HLA-DQA2")) %>% 
  filter(!(gene == 3802 & symbol == "KIR2DS1")) %>% 
  filter(!(gene == 414243 & symbol == "LINC00856")) %>% 
  filter(!(gene == 100132596 & symbol == "XG"))

# genes in analysis
cg = read_tsv("~/projects/age-sex-prediction/data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", col_names=F) %>% pull(1)
# fix case of gene symbols
u <- sc %>% filter(toupper(geneSymbolmore1) %in% gene_sym$symbol) %>% pull(geneSymbolmore1) %>% unique
l <- sc %>% filter(geneSymbolmore1 %in% gene_sym$symbol) %>% pull(geneSymbolmore1) %>% unique
lgenes <- setdiff(l, u)

sc <- sc %>% 
  mutate(geneSymbolmore1 = ifelse(geneSymbolmore1 %in% lgenes, geneSymbolmore1, toupper(geneSymbolmore1))) %>% 
  rename(symbol = geneSymbolmore1)

sc <- left_join(sc, gene_sym, by = "symbol")

sc %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sctype/sctype_gene_info.tsv")

sc %>% 
  unite(group, tissueType, cellName, sep = ": ") %>% 
  select(gene, group) %>% 
  rename(entrez = gene) %>% 
  filter(!is.na(entrez)) %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sctype/sctype_gene-group_file.tsv")

print("script finished")

