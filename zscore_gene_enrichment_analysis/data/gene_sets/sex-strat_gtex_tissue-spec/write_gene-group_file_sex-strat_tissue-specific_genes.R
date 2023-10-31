tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

# Function requires a vector with expression of one gene in different tissues.
fTau <- function(x){
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      # print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    # print("No data for this gene available.")
  } 
  return(res)
}

# Function requires a vector with expression of one gene in different tissues.
robust_z <- function(x){
  med <- median(x)
  rz <- (x - med)/(median(abs(x-med)) * 1.4826)
  return(rz)
}

gtex <- read_delim("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sex-strat_gtex_tissue-spec/GTEX_median_tpm_for_each_tissue_by_sex.tsv",
                   delim = "\t", col_names = T)

gtex <- gtex %>% 
  mutate(median_tpm = asinh(median_tpm))

out <- tibble()

# female tau/robust z
fgtex <- gtex %>% 
  filter(sex == "female")

for (gene in unique(fgtex$Entrez)){
  gene_df <- fgtex %>% 
    filter(Entrez == gene)
  
  gene_tau <- fTau(pull(gene_df, median_tpm))
  
  gene_df <- gene_df %>% 
    mutate(robust_zscore = robust_z(median_tpm)) %>% 
    mutate(tau_score = gene_tau)
  
  out <- bind_rows(out, gene_df)
}

# male tau/robust z
mgtex <- gtex %>% 
  filter(sex == "male")

for (gene in unique(mgtex$Entrez)){
  gene_df <- mgtex %>% 
    filter(Entrez == gene)
  
  gene_tau <- fTau(pull(gene_df, median_tpm))
  
  gene_df <- gene_df %>% 
    mutate(robust_zscore = robust_z(median_tpm)) %>% 
    mutate(tau_score = gene_tau)
  
  out <- bind_rows(out, gene_df)
}

# write tau/robust z scores
out %>% 
  write_tsv("./GTEX_tau_robust-z_scores_tissue_each_sex.tsv")

### write tissue specific genes for each sex ###
out %>% 
  filter(tau_score > 0.8) %>% 
  filter(robust_zscore >= 2) %>% 
  mutate(tissue = paste(sex, tissue, sep = "__")) %>% 
  select(Entrez, tissue) %>% 
  write_tsv("./sex-strat_gtex_tissue-spec_gene-group_file.tsv")

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

