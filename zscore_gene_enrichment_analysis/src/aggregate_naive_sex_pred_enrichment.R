tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)
# args[1] = gene_set, one of:
#                       disgenetC
#                       disgenetF
#                       gobp
#                       gtex
#                       guo
#                       gwas
#                       mp
#                       sagd
gene_set <- args[1]

results_dir <- "../results/naive_sex_prediction_ranks/"
rnaseq_results_dir <- paste0(results_dir, "rnaseq/", gene_set)
microarray_results_dir <- paste0(results_dir, "microarray/", gene_set)

# results files cols: group | zscore
# list files
rnaseq_files <- list.files(rnaseq_results_dir, full.names = T)
# read into list
rnaseq_results <- lapply(rnaseq_files, read_tsv)
# name list for id col
names(rnaseq_results) <- basename(rnaseq_files)
# bind list of dfs into one df
rnaseq_results <- bind_rows(rnaseq_results, .id = "file")
# create age group column
rnaseq_results <- rnaseq_results %>% 
  separate(file, into = c("age_group", "file"), sep = "___") %>% 
  select(age_group, group, zscore)
colnames(rnaseq_results) <- c("age_group", "term", "zscore")

# list files
microarray_files <- list.files(microarray_results_dir, full.names = T)
# read into list
microarray_results <- lapply(microarray_files, read_tsv)
# name list for id col
names(microarray_results) <- basename(microarray_files)
# bind list of dfs into one df
microarray_results <- bind_rows(microarray_results, .id = "file")
# create age group column
microarray_results <- microarray_results %>% 
  separate(file, into = c("age_group", "file"), sep = "___") %>% 
  select(age_group, group, zscore)
colnames(microarray_results) <- c("age_group", "term", "zscore")

# combine results
results <- bind_rows(rnaseq_results, microarray_results)
results <- results %>%  
  group_by(age_group, term) %>% 
  mutate(combined_zscore = (sum(zscore) / sqrt(2))) %>% 
  select(age_group, term, combined_zscore) %>% 
  distinct() %>% 
  ungroup()

results %>% write_tsv(paste0(results_dir, gene_set, "_enrichment.tsv"))

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
