library(tidyverse)
# this script puts together all output files for analysis

auroc_files <- list.files("../../results/age_prediction_by_sex/microarray/output_files/",
                          full.names = T, pattern = "*_auroc_*")
sample_prob_files <- list.files("../../results/age_prediction_by_sex/microarray/output_files/",
                                full.names = T, pattern = "*_sample_probabilities_*")
model_weights_files <- list.files("../../results/age_prediction_by_sex/microarray/output_files/",
                                  full.names = T, pattern = "*_model_weights_*")

auroc <- lapply(auroc_files, read_delim, delim = "\t", col_names = T)
names(auroc) <- basename(auroc_files)
auroc <- bind_rows(auroc, .id = "parameters")

auroc <- auroc %>% 
  mutate(penalty = ifelse(grepl("l2_LR", parameters), "l2", "elasticnet")) %>% 
  mutate(std_scld = ifelse(grepl("without_std_scaling", parameters), F, T)) %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(group_type = ifelse(grepl("coarse", parameters), "coarse", "fine")) %>% 
  mutate(age_group = ifelse(grepl("fetus", parameters), "fetus",
                            ifelse(grepl("infant", parameters), "infant",
                                   ifelse(grepl("young_child", parameters), "young_child",
                                          ifelse(grepl("ale_child_results", parameters), "child",
                                                 ifelse(grepl("adolescent", parameters), "adolescent",
                                                        ifelse(grepl("young_adult", parameters), "young_adult",
                                                               ifelse(grepl("ale_adult_results", parameters), "adult",
                                                                      ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                             ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                    ifelse(grepl("old_adult", parameters), "old_adult", "elderly")))))))))))

sample_probs <- lapply(sample_prob_files, read_delim, delim = "\t", col_names = T)
names(sample_probs) <- basename(sample_prob_files)
sample_probs <- lapply(sample_probs, mutate, pos_age_group = as.character(pos_age_group))
sample_probs <- bind_rows(sample_probs, .id = "parameters")

sample_probs <- sample_probs %>% 
  mutate(penalty = ifelse(grepl("l2_LR", parameters), "l2", "elasticnet")) %>% 
  mutate(std_scld = ifelse(grepl("without_std_scaling", parameters), F, T)) %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(group_type = ifelse(grepl("coarse", parameters), "coarse", "fine")) %>% 
  mutate(age_group = ifelse(grepl("fetus", parameters), "fetus",
                            ifelse(grepl("infant", parameters), "infant",
                                   ifelse(grepl("young_child", parameters), "young_child",
                                          ifelse(grepl("ale_child_results", parameters), "child",
                                                 ifelse(grepl("adolescent", parameters), "adolescent",
                                                        ifelse(grepl("young_adult", parameters), "young_adult",
                                                               ifelse(grepl("ale_adult_results", parameters), "adult",
                                                                      ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                             ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                    ifelse(grepl("old_adult", parameters), "old_adult", "elderly"))))))))))) %>% 
  filter(penalty == "elasticnet") %>% 
  filter(std_scld == FALSE)


model_weights <- lapply(model_weights_files, read_delim, delim = "\t", col_names = T)
names(model_weights) <- basename(model_weights_files)
model_weights <- bind_rows(model_weights, .id = "parameters")

model_weights <- model_weights %>% 
  mutate(penalty = ifelse(grepl("l2_LR", parameters), "l2", "elasticnet")) %>% 
  mutate(std_scld = ifelse(grepl("without_std_scaling", parameters), F, T)) %>% 
  mutate(sex = ifelse(grepl("female", parameters), "female", "male")) %>% 
  mutate(group_type = ifelse(grepl("coarse", parameters), "coarse", "fine")) %>% 
  mutate(age_group = ifelse(grepl("fetus", parameters), "fetus",
                            ifelse(grepl("infant", parameters), "infant",
                                   ifelse(grepl("young_child", parameters), "young_child",
                                          ifelse(grepl("ale_child_weights", parameters), "child",
                                                 ifelse(grepl("adolescent", parameters), "adolescent",
                                                        ifelse(grepl("young_adult", parameters), "young_adult",
                                                               ifelse(grepl("ale_adult_weights", parameters), "adult",
                                                                      ifelse(grepl("middle_adult", parameters), "middle_adult",
                                                                             ifelse(grepl("older_adult", parameters), "older_adult",
                                                                                    ifelse(grepl("old_adult", parameters), "old_adult", "elderly"))))))))))) %>% 
  filter(penalty == "elasticnet") %>% 
  filter(std_scld == FALSE)

auroc %>% 
  write_delim("../../results/age_prediction_by_sex/microarray/auroc_results.tsv",
              delim = "\t", col_names = T)

sample_probs %>% 
  write_delim("../../results/age_prediction_by_sex/microarray/sample_probabilites_from_model_with_chosen_parameters.tsv",
              delim = "\t", col_names = T)

model_weights %>% 
  write_delim("../../results/age_prediction_by_sex/microarray/model_weights_from_model_with_chosen_parameters.tsv",
              delim = "\t", col_names = T)