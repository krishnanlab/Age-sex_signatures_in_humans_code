tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

file_path <- "~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sex-strat_gtex_tissue-spec/ind_sex-tissue_outputs"

file_list <- list.files(file_path, full.names = T)

# read
files <- lapply(file_list, read_tsv)

data <- bind_rows(files)

# remove ensg column, avg any multimapped genes
data <- data %>% 
  filter(!is.na(Entrez)) %>%
  select(-ENSG) %>% 
  distinct() %>% 
  group_by(Entrez, sex, tissue) %>% 
  summarise(median_tpm = mean(median_tpm)) %>% 
  ungroup()

# write file
data %>% 
  write_tsv("./GTEX_median_tpm_for_each_tissue_by_sex.tsv")

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
