library(tidyverse)

sex_dir <- "~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/naive_sex_prediction_signed_ranks/"
age_dir <- "~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/elasticnet_LR_model_weights/"

ar_tibble <- tibble(age_group = c("fetus","infant","young_child", "child","adolescent", "young_adult",
                                  "adult", "middle_adult", "older_adult", "old_adult", "elderly"),
                    age_range = c("< 0", "[0-2]", "(2-8]", "(8-12]", "(12-20]", "(20-35]",
                                  "(35-45]", "(45-60]", "(60-70]", "(70-80]", "> 80"))

age_ranges <- c("< 0", "[0-2]", "(2-8]", "(8-12]", "(12-20]", "(20-35]",
                "(35-45]", "(45-60]", "(60-70]", "(70-80]", "> 80")

# put scores together for shiny app
# separately for RNAseq/microarray

# sex
file_list <- list.files(paste0(sex_dir, "rnaseq"), full.names = T)

out <- tibble()

for(file in file_list){
  fname <- gsub("_signed_gene_cut_data_all_chromosomes.tsv", "", basename(file))
  dat <- read_tsv(file) %>% 
    mutate(age_group = fname) %>% 
    mutate(technology = "RNAseq")
  out <- bind_rows(out, dat)
}

file_list <- list.files(paste0(sex_dir, "microarray"), full.names = T)

for(file in file_list){
  fname <- gsub("_signed_gene_cut_data_all_chromosomes.tsv", "", basename(file))
  dat <- read_tsv(file) %>% 
    mutate(age_group = fname) %>% 
    mutate(technology = "microarray")
  out <- bind_rows(out, dat)
}

out <- left_join(out, ar_tibble, by = "age_group")

out %>% 
  rename(score = signed_balanced_accuracy) %>% 
  write_tsv("~/projects/Age-sex_signatures_in_humans_app/data/sex_bias_scores.tsv")


# age
file_list <- list.files(paste0(age_dir, "rnaseq"), pattern = "*fine*", full.names = T)

out <- tibble()

for(file in file_list){
  fname <- gsub("_RNAseq_model_weights.tsv", "", basename(file))
  fname <- gsub("-fine-1", "", fname)
  fname <- gsub("-fine-2", "", fname)
  fname <- gsub("-fine-3", "", fname)
  dat <- read_tsv(file) 
  colnames(dat) <- c("gene", "score")
  dat <- dat %>% 
    mutate(tmp = fname) %>% 
    separate_wider_delim(tmp, names = c("sex", "age_group"), delim = "-") %>% 
    mutate(technology = "RNAseq") 
  out <- bind_rows(out, dat)
}

file_list <- list.files(paste0(age_dir, "microarray"), pattern = "*fine*", full.names = T)

for(file in file_list){
  fname <- gsub("_microarray_model_weights.tsv", "", basename(file))
  fname <- gsub("-fine-1", "", fname)
  fname <- gsub("-fine-2", "", fname)
  fname <- gsub("-fine-3", "", fname)
  dat <- read_tsv(file) 
  colnames(dat) <- c("gene", "score")
  dat <- dat %>% 
    mutate(tmp = fname) %>% 
    separate_wider_delim(tmp, names = c("sex", "age_group"), delim = "-") %>%  
    mutate(technology = "microarray")
  out <- bind_rows(out, dat)
}

out <- out %>% 
  group_by(gene, sex, age_group, technology) %>% 
  dplyr::summarise(score = median(score))

out <- left_join(out, ar_tibble, by = "age_group")

out %>% 
  write_tsv("~/projects/Age-sex_signatures_in_humans_app/data/age_bias_scores.tsv")


