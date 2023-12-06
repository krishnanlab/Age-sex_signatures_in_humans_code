library(tidyverse)
# this script reads the all chromosome naive sex pred results in and 
# makes a "signed_balanced_accuracy" column which puts a negative sign
# on the balanced accuracy of genes which have higher expression in males 
# then writes a new file for each age group

rf_path <- "../../zscore_gene_enrichment_analysis/data/naive_sex_prediction_signed_ranks/rnaseq/"
mf_path <- "../../zscore_gene_enrichment_analysis/data/naive_sex_prediction_signed_ranks/microarray/"

rf <- read_tsv("../../results/naive_sex_prediction/rnaseq/final_gene_cut_data_all_chromosomes.tsv")
mf <- read_tsv("../../results/naive_sex_prediction/microarray/final_gene_cut_data_all_chromosomes.tsv")

# fix genes with no cutoff
rf <- rf %>% 
  mutate(balanced_accuracy = ifelse(cutoff == "none", 0.5, balanced_accuracy)) %>% 
  mutate(higher_sex = ifelse(cutoff == "none", "none", higher_sex)) %>% 
  mutate(cutoff = ifelse(cutoff == "none", 0, cutoff)) %>% 
  mutate(cutoff = as.numeric(cutoff))

common_genes <- read_tsv("../../data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", col_names = F) %>% 
  pull(1)

rf <- rf %>% 
  filter(gene %in% common_genes) %>% 
  mutate(fhigh_balanced_accuracy = ifelse(higher_sex == "male", 1-balanced_accuracy, balanced_accuracy))

mf <- mf %>% 
  filter(gene %in% common_genes) %>% 
  mutate(fhigh_balanced_accuracy = ifelse(higher_sex == "male", 1-balanced_accuracy, balanced_accuracy)) 

fine_ags <- c("fetus", "infant", "young_child", "child", "adolescent", "young_adult", 
              "adult", "middle_adult", "older_adult", "old_adult", "elderly")

no_var_genes <- setdiff(common_genes, unique(rf$gene))

rfaddition <-  tibble(age_group = rep(fine_ags, each = 164), gene = rep(no_var_genes, times = 11), cutoff = 0, 
                      balanced_accuracy = 0.5, higher_sex = 'none', 
                      fhigh_balanced_accuracy = 0.5)

rf <- bind_rows(rf, rfaddition)

rf <- rf %>% 
  mutate(signed_balanced_accuracy = ((2 * fhigh_balanced_accuracy) - 1))
mf <- mf %>% 
  mutate(signed_balanced_accuracy = ((2 * fhigh_balanced_accuracy) - 1))

for (ag in fine_ags){
  rf %>% 
    filter(age_group == ag) %>% 
    select(gene, signed_balanced_accuracy) %>% 
    write_delim(paste0(rf_path, ag, "_signed_gene_cut_data_all_chromosomes.tsv"),
                delim = "\t", col_names = T)
}

for (ag in fine_ags){
  mf %>% 
    filter(age_group == ag) %>% 
    select(gene, signed_balanced_accuracy) %>% 
    write_delim(paste0(mf_path, ag, "_signed_gene_cut_data_all_chromosomes.tsv"),
                delim = "\t", col_names = T)
}




