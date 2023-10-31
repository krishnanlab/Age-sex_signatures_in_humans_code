library(tidyverse)
# this script puts together all output files for analysis

auroc_files <- list.files("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/output_files/",
                          full.names = T, pattern = "*_auroc_*")
auprc_files <- list.files("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/output_files/",
                          full.names = T, pattern = "*_auprc_*")
sample_prob_files <- list.files("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/output_files/",
                                full.names = T, pattern = "*_sample_probabilities_*")
model_weights_files <- list.files("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/output_files/",
                                  full.names = T, pattern = "*_model_weights_*")

# read files, bind into one df
auroc <- lapply(auroc_files, read_delim, delim = "\t", col_names = T)
names(auroc) <- basename(auroc_files)
auroc <- bind_rows(auroc, .id = "parameters")
# pull sex/age out of file names
auroc <- auroc %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(model_age_group = ifelse(grepl("fetus", parameters), "fetus",
                                  ifelse(grepl("infant", parameters), "infant",
                                         ifelse(grepl("young_child", parameters), "young_child",
                                                ifelse(grepl("ale_child_results", parameters), "child",
                                                       ifelse(grepl("adolescent", parameters), "adolescent",
                                                              ifelse(grepl("young_adult", parameters), "young_adult",
                                                                     ifelse(grepl("ale_adult_results", parameters), "adult",
                                                                            ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                                   ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                          ifelse(grepl("old_adult", parameters), "old_adult", "elderly")))))))))))
# read files, bind into one df
auprc <- lapply(auprc_files, read_delim, delim = "\t", col_names = T)
names(auprc) <- basename(auprc_files)
auprc <- bind_rows(auprc, .id = "parameters")
# pull sex/age out of file names
auprc <- auprc %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(model_age_group = ifelse(grepl("fetus", parameters), "fetus",
                                  ifelse(grepl("infant", parameters), "infant",
                                         ifelse(grepl("young_child", parameters), "young_child",
                                                ifelse(grepl("ale_child_results", parameters), "child",
                                                       ifelse(grepl("adolescent", parameters), "adolescent",
                                                              ifelse(grepl("young_adult", parameters), "young_adult",
                                                                     ifelse(grepl("ale_adult_results", parameters), "adult",
                                                                            ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                                   ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                          ifelse(grepl("old_adult", parameters), "old_adult", "elderly")))))))))))

# read files, bind into one df
sample_probs <- lapply(sample_prob_files, read_delim, delim = "\t", col_names = T)
names(sample_probs) <- basename(sample_prob_files)
sample_probs <- bind_rows(sample_probs, .id = "parameters")

sample_probs <- sample_probs %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(model_age_group = ifelse(grepl("fetus", parameters), "fetus",
                            ifelse(grepl("infant", parameters), "infant",
                                   ifelse(grepl("young_child", parameters), "young_child",
                                          ifelse(grepl("ale_child_results", parameters), "child",
                                                 ifelse(grepl("adolescent", parameters), "adolescent",
                                                        ifelse(grepl("young_adult", parameters), "young_adult",
                                                               ifelse(grepl("ale_adult_results", parameters), "adult",
                                                                      ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                             ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                    ifelse(grepl("old_adult", parameters), "old_adult", "elderly")))))))))))


# read files, bind into one df
model_weights <- lapply(model_weights_files, read_delim, delim = "\t", col_names = T)
names(model_weights) <- basename(model_weights_files)
model_weights <- bind_rows(model_weights, .id = "parameters")

model_weights <- model_weights %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(model_age_group = ifelse(grepl("fetus", parameters), "fetus",
                            ifelse(grepl("infant", parameters), "infant",
                                   ifelse(grepl("young_child", parameters), "young_child",
                                          ifelse(grepl("ale_child_weights", parameters), "child",
                                                 ifelse(grepl("adolescent", parameters), "adolescent",
                                                        ifelse(grepl("young_adult", parameters), "young_adult",
                                                               ifelse(grepl("ale_adult_weights", parameters), "adult",
                                                                      ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                             ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                    ifelse(grepl("old_adult", parameters), "old_adult", "elderly")))))))))))

# write condensed results
auroc %>% 
  write_delim("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/auroc_results_from_rnaseq_blood_models.tsv",
              delim = "\t", col_names = T)

auprc %>% 
  write_delim("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/auprc_results_from_rnaseq_blood_models.tsv",
              delim = "\t", col_names = T)

sample_probs %>% 
  write_delim("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/sample_probabilites_from_rnaseq_blood_models.tsv",
              delim = "\t", col_names = T)

model_weights %>% 
  write_delim("~/projects/age-sex-prediction/blood-only_analysis/results/age_prediction_by_sex/rnaseq/model_weights_from_rnaseq_blood_models.tsv",
              delim = "\t", col_names = T)


print("done")
