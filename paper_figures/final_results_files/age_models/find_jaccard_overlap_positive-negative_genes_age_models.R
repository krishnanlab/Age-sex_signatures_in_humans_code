library(tidyverse)

# read data
rmodel_weights <- read_delim("~/projects/age-sex-prediction/results/age_prediction_by_sex/rnaseq/model_weights_from_model_with_chosen_parameters.tsv", 
                             delim="\t", col_names=T) %>% 
  filter(group_type == "fine") %>%
  select(-std_scld, -penalty, -parameters, -asinh, -n_positives, -group_type) %>% 
  select(sex, age_group, everything())

mmodel_weights <- read_delim("~/projects/age-sex-prediction/results/age_prediction_by_sex/microarray/model_weights_from_model_with_chosen_parameters.tsv", 
                             delim="\t", col_names=T) %>% 
  filter(group_type == "fine") %>%
  select(-std_scld, -penalty, -parameters, -n_positives, -group_type) %>% 
  select(sex, age_group, everything())

mmodel_weights <- mmodel_weights[names(rmodel_weights)]

# add fold to keep unique
rmodel_weights <- rmodel_weights %>% 
  mutate(age_group = paste(age_group, c(1,2,3), sep = "_"))
# add fold to keep unique
mmodel_weights <- mmodel_weights %>% 
  mutate(age_group = paste(age_group, c(1,2,3), sep = "_"))
# make tidy
rmodel_weights <- rmodel_weights %>% 
  pivot_longer(3:18480, names_to = "gene", values_to = "weight")

mmodel_weights <- mmodel_weights %>% 
  pivot_longer(3:18480, names_to = "gene", values_to = "weight")

# add gene grouping col/data col
rmodel_weights <- rmodel_weights %>% 
  mutate(gene_group = ifelse(weight > 0, "positive",
                             ifelse(weight < 0, "negative", "zero"))) %>% 
  mutate(data = "RNAseq")

mmodel_weights <- mmodel_weights %>% 
  mutate(gene_group = ifelse(weight > 0, "positive",
                             ifelse(weight < 0, "negative", "zero"))) %>% 
  mutate(data = "microarray")

model_weights <- bind_rows(rmodel_weights, mmodel_weights)

ab_positive_intersect <- function(df){
  df <- df %>% 
    unite(group, data, sex, age_group, sep = "_")
  
  out <- tibble(group_one = "group",
                group_two = "group",
                intersection = 1,
                union = 2)
  
  for (a in unique(df$group)){
    for (b in unique(df$group)){
      a_genes <- df %>% 
        filter(gene_group == "positive") %>% 
        filter(group == !!a) %>% 
        pull(gene) %>% 
        unique()
      b_genes <- df %>% 
        filter(gene_group == "positive") %>% 
        filter(group == !!b) %>% 
        pull(gene) %>% 
        unique()
      i <- length(intersect(a_genes, b_genes))
      u <- length(union(a_genes, b_genes))
      tmp <- tibble(group_one = a,
                    group_two = b,
                    intersection = i,
                    union = u)
      out <- bind_rows(out, tmp)
    }
  }
  
  out <- out[-1,]
  out <- out %>% 
    mutate(jaccard = intersection / union)
  return(out)  
}

ab_negative_intersect <- function(df){
  df <- df %>% 
    unite(group, data, sex, age_group, sep = "_")
  
  out <- tibble(group_one = "group",
                group_two = "group",
                intersection = 1,
                union = 2)
  
  for (a in unique(df$group)){
    for (b in unique(df$group)){
      a_genes <- df %>% 
        filter(gene_group == "negative") %>% 
        filter(group == !!a) %>% 
        pull(gene) %>% 
        unique()
      b_genes <- df %>% 
        filter(gene_group == "negative") %>% 
        filter(group == !!b) %>% 
        pull(gene) %>% 
        unique()
      i <- length(intersect(a_genes, b_genes))
      u <- length(union(a_genes, b_genes))
      tmp <- tibble(group_one = a,
                    group_two = b,
                    intersection = i,
                    union = u)
      out <- bind_rows(out, tmp)
    }
  }
  
  out <- out[-1,]
  out <- out %>% 
    mutate(jaccard = intersection / union)
  return(out)  
}

pos_int <- ab_positive_intersect(model_weights)

neg_int <- ab_negative_intersect(model_weights)

pos_int %>%
  write_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/age_models/jaccard_overlap_of_positive_genes_age_models.tsv")

neg_int %>%
  write_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/age_models/jaccard_overlap_of_negative_genes_age_models.tsv")

