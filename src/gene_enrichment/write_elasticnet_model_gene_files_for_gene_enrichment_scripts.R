library(tidyverse)

tic <- as.integer(as.POSIXct(Sys.time()))

######## RNAseq ########
# n_positives is the # of positives in the test fold
model_weights <- read_delim("~/projects/age-sex-prediction/results/age_prediction_by_sex/rnaseq/model_weights_from_model_with_chosen_parameters.tsv", 
                            delim="\t", col_names=T) %>% 
  mutate(model_group = paste(sex, age_group, group_type, c("1", "2", "3"), sep = "-")) %>% 
  select(-std_scld, -penalty, -parameters, -asinh, -n_positives, -age_group, -sex, -group_type) %>% 
  select(model_group, everything())
  
# transpose so easier to work with
new_cols <- model_weights %>% pull(model_group)
model_weights <- model_weights %>% 
  select(-model_group) %>% 
  as.matrix()

model_weights <- t(model_weights)
genes <- rownames(model_weights)
rownames(model_weights) <- NULL

# zscore
zmodel_weights <- scale(model_weights)
colnames(zmodel_weights) <- new_cols
zmodel_weights <- as_tibble(zmodel_weights) %>% 
  mutate(gene = genes) %>% 
  select(gene, everything())

# add names back to not z scored wts
colnames(model_weights) <- new_cols
model_weights <- as_tibble(model_weights) %>% 
  mutate(gene = genes) %>% 
  select(gene, everything())

# get model groups (but not 'gene')
groups <- colnames(zmodel_weights)
groups <- groups[-1]

output <- tibble(gene = "1", group = "group")
for (group in groups) {
  df <- zmodel_weights %>% 
    select(gene, all_of(group))
 neg_genes <-  df %>% 
    filter(!!as.name(group) < -2) %>% 
    pull(gene)
  pos_genes <- df %>% 
    filter(!!as.name(group) > 2) %>% 
    pull(gene)
  pos_group_df <- tibble(gene = pos_genes, group = paste0("pos_", group))
  neg_group_df <- tibble(gene = neg_genes, group = paste0("neg_", group))
  group_df <- bind_rows(pos_group_df, neg_group_df)
  output <- bind_rows(output, group_df)
}
output <- output[-1,]

output %>% 
  write_delim("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/RNAseq_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv", 
              delim = "\t", col_names = T)

for (group in groups) {
  df <- model_weights %>% 
    select(gene, all_of(group)) %>% 
    write_delim(paste0("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/",
                       group, "_RNAseq_model_weights.tsv"),
                delim = "\t", col_names = T)
}

######## Microarray ########
# n_positives is the # of positives in the test fold
model_weights <- read_delim("~/projects/age-sex-prediction/results/age_prediction_by_sex/microarray/model_weights_from_model_with_chosen_parameters.tsv", 
                            delim="\t", col_names=T) %>% 
  mutate(model_group = paste(sex, age_group, group_type, c("1", "2", "3"), sep = "-")) %>% 
  select(-std_scld, -penalty, -parameters, -n_positives, -age_group, -sex, -group_type) %>% 
  select(model_group, everything())

# transpose so easier to work with
new_cols <- model_weights %>% pull(model_group)
model_weights <- model_weights %>% 
  select(-model_group) %>% 
  as.matrix()

model_weights <- t(model_weights)
genes <- rownames(model_weights)
rownames(model_weights) <- NULL

# zscore
zmodel_weights <- scale(model_weights)
colnames(zmodel_weights) <- new_cols
zmodel_weights <- as_tibble(zmodel_weights) %>% 
  mutate(gene = genes) %>% 
  select(gene, everything())

# add names back to not z scored wts
colnames(model_weights) <- new_cols
model_weights <- as_tibble(model_weights) %>% 
  mutate(gene = genes) %>% 
  select(gene, everything())

# get model groups (but not 'gene')
groups <- colnames(zmodel_weights)
groups <- groups[-1]

output <- tibble(gene = "1", group = "group")
for (group in groups) {
  df <- zmodel_weights %>% 
    select(gene, all_of(group))
  neg_genes <-  df %>% 
    filter(!!as.name(group) < -2) %>% 
    pull(gene)
  pos_genes <- df %>% 
    filter(!!as.name(group) > 2) %>% 
    pull(gene)
  pos_group_df <- tibble(gene = pos_genes, group = paste0("pos_", group))
  neg_group_df <- tibble(gene = neg_genes, group = paste0("neg_", group))
  group_df <- bind_rows(pos_group_df, neg_group_df)
  output <- bind_rows(output, group_df)
}
output <- output[-1,]

output %>% 
  write_delim("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/microarray_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv", 
              delim = "\t", col_names = T)

for (group in groups) {
  df <- model_weights %>% 
    select(gene, all_of(group)) %>% 
    write_delim(paste0("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/",
                       group, "_microarray_model_weights.tsv"),
                delim = "\t", col_names = T)
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
