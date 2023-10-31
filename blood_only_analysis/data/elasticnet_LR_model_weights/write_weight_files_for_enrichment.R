library(tidyverse)

############ RNAseq ############
weights <- read_tsv("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/model_weights_from_rnaseq_blood_models.tsv")
# get rid of unnecessary parts of file name (keep sex/age group)
weights$parameters <- gsub("3FCV_elasticnet_LR_model_weights_asinh_transformed_no_std_scaling_", "", weights$parameters)
weights$parameters <- gsub("_weights_with_fine_ags.tsv", "", weights$parameters)
# add fold
weights$parameters <- paste(weights$parameters, rep(c("1","2","3"), times = 18), sep = "_")

# models
models <- weights %>% pull(parameters)
# get rid of nonweight cols
weights <- weights %>% 
  select(-parameters) %>% 
  select(-n_positives) %>% 
  select(-sex) %>% 
  select(-model_age_group)
# genes, get rid of "parameters"
genes <- colnames(weights)

# transpose
colnames(weights) <- NULL
weights <- as.matrix(weights)
weights <- t(weights)
colnames(weights) <- models
weights <- as_tibble(weights)
weights$gene <- genes

# write files
out_path <- "~/projects/age-sex-prediction/blood-only_analysis/data/elasticnet_LR_model_weights/rnaseq/"
for (model in models){
  out <- weights %>% 
    select(gene, !!model)
  out %>% 
    write_tsv(paste0(out_path, model, "_weights.tsv"))
}

############ Microarray ############
weights <- read_tsv("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/microarray/model_weights_from_microarray_blood_models.tsv")
# get rid of unnecessary parts of file name (keep sex/age group)
weights$parameters <- gsub("3FCV_elasticnet_LR_model_weights_no_std_scaling_", "", weights$parameters)
weights$parameters <- gsub("_weights_with_fine_ags.tsv", "", weights$parameters)
# add fold
weights$parameters <- paste(weights$parameters, rep(c("1","2","3"), times = 18), sep = "_")

# models
models <- weights %>% pull(parameters)
# get rid of nonweight cols
weights <- weights %>% 
  select(-parameters) %>% 
  select(-n_positives) %>% 
  select(-sex) %>% 
  select(-model_age_group)
# genes, get rid of "parameters"
genes <- colnames(weights)

# transpose
colnames(weights) <- NULL
weights <- as.matrix(weights)
weights <- t(weights)
colnames(weights) <- models
weights <- as_tibble(weights)
weights$gene <- genes

# write files
out_path <- "~/projects/age-sex-prediction/blood-only_analysis/data/elasticnet_LR_model_weights/microarray/"
for (model in models){
  out <- weights %>% 
    select(gene, !!model)
  out %>% 
    write_tsv(paste0(out_path, model, "_weights.tsv"))
}


print("done")
