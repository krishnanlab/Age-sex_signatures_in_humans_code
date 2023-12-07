tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

enrichment_dir <- "../../results/naive_sex_prediction_ranks/"
enrichment_files <- list.files(enrichment_dir, pattern = "tsv", full.names = T)
enrichment_files <- enrichment_files[!grepl("disgenet", enrichment_files)]
output_dir <- "../../data/orsum_blood_naive_sex_prediction_term_significance_files/"

# write terms in order of significance for each age group and term set
for (ef in enrichment_files){
  res <- read_tsv(ef)
  ags <- unique(res$age_group)
  term_set <- basename(ef)
  term_set <- gsub("_enrichment_blood_samples.tsv", "", term_set)
  for (ag in ags){
    res %>% 
      filter(age_group == ag) %>% 
      arrange(desc(abs(combined_zscore))) %>% 
      select(term) %>% 
      write_tsv(paste0(output_dir, term_set, "_", ag, "_terms_by_significance.txt"),
                col_names = F)
  }
}

enrichment_dir <- "../../results/elasticnet_LR_model_weights/"
enrichment_files <- list.files(enrichment_dir, pattern = "tsv", full.names = T)
enrichment_files <- enrichment_files[!grepl("disgenet", enrichment_files)]
output_dir <- "../../data/orsum_blood_age_prediction_significance_files/"

# write terms in order of significance for each sex-age group and term set
for (ef in enrichment_files){
  res <- read_tsv(ef)
  ags <- unique(res$age_group)
  term_set <- basename(ef)
  term_set <- gsub("_enrichment_blood_samples.tsv", "", term_set)
  for (ag in ags){
    for (s in c("female", "male")){
      res %>% 
        filter(age_group == ag) %>% 
        filter(sex == s) %>% 
        arrange(desc(abs(combined_zscore))) %>% 
        select(term) %>% 
        write_tsv(paste0(output_dir, term_set, "_", s, "_", ag, "_terms_by_significance.txt"),
                  col_names = F)
    }
  }
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
