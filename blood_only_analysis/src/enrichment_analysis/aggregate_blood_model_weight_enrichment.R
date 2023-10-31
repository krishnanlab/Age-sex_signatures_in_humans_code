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

results_dir <- "~/projects/age-sex-prediction/blood-only_analysis/results/elasticnet_LR_model_weights/"
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
# create sex/age group column
rnaseq_results <- rnaseq_results %>% 
  separate(file, into = c("age_group", "file"), sep = "___") %>% 
  mutate(sex = ifelse(grepl("female", age_group), "female", "male"))
rnaseq_results$age_group <- gsub("female_", "", rnaseq_results$age_group)
rnaseq_results$age_group <- gsub("male_", "", rnaseq_results$age_group)
rnaseq_results$age_group <- gsub("_1", "", rnaseq_results$age_group)
rnaseq_results$age_group <- gsub("_2", "", rnaseq_results$age_group)
rnaseq_results$age_group <- gsub("_3", "", rnaseq_results$age_group)
rnaseq_results <- rnaseq_results %>% 
  select(sex, age_group, group, zscore)
colnames(rnaseq_results) <- c("sex", "age_group", "term", "zscore")
rnaseq_results <- rnaseq_results %>% 
  group_by(sex, age_group, term) %>% 
  mutate(avg_zscore = mean(zscore)) %>% 
  select(sex, age_group, term, avg_zscore) %>% 
  distinct() %>% 
  ungroup()

# list files
microarray_files <- list.files(microarray_results_dir, full.names = T)
# read into list
microarray_results <- lapply(microarray_files, read_tsv)
# name list for id col
names(microarray_results) <- basename(microarray_files)
# bind list of dfs into one df
microarray_results <- bind_rows(microarray_results, .id = "file")
# create sex/age group column
microarray_results <- microarray_results %>% 
  separate(file, into = c("age_group", "file"), sep = "___") %>% 
  mutate(sex = ifelse(grepl("female", age_group), "female", "male"))
microarray_results$age_group <- gsub("female_", "", microarray_results$age_group)
microarray_results$age_group <- gsub("male_", "", microarray_results$age_group)
microarray_results$age_group <- gsub("_1", "", microarray_results$age_group)
microarray_results$age_group <- gsub("_2", "", microarray_results$age_group)
microarray_results$age_group <- gsub("_3", "", microarray_results$age_group)
microarray_results <- microarray_results %>% 
  select(sex, age_group, group, zscore)
colnames(microarray_results) <- c("sex", "age_group", "term", "zscore")
microarray_results <- microarray_results %>% 
  group_by(sex, age_group, term) %>% 
  mutate(avg_zscore = mean(zscore)) %>% 
  select(sex, age_group, term, avg_zscore) %>% 
  distinct() %>% 
  ungroup()

# combine results
results <- bind_rows(rnaseq_results, microarray_results)
results <- results %>%  
  group_by(sex, age_group, term) %>% 
  mutate(combined_zscore = (sum(avg_zscore) / sqrt(2))) %>% 
  select(sex, age_group, term, combined_zscore) %>% 
  distinct() %>% 
  ungroup()

results %>% write_tsv(paste0(results_dir, gene_set, "_enrichment_blood_samples.tsv"))

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
