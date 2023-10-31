library(tidyverse)
args <- commandArgs(TRUE)

tic <- as.integer(as.POSIXct(Sys.time()))

out_path <- "~/projects/age-sex-prediction/results/age_prediction_by_sex/gene_enrichment_analysis/output_files/"
# data_path <- "~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/"

source("findEnrichedTerms.R")

num_common_genes = read_delim("~/projects/age-sex-prediction/data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", 
                          delim = "\t", col_names = F) %>% 
  pull(1) %>% 
  length()

set_one <- read_delim(args[1], delim = "\t", col_names = T)
set_two <- read_delim(args[2], delim = "\t", col_names = T)

output <- overlapSets(table1 = set_one,
                      table2 = set_two,
                      background = num_common_genes)

set_one_name <- basename(args[1])
set_one_name <- gsub("_gene-group_file.tsv", "", set_one_name)
set_one_name <- gsub("_genes_over_2sd_from_mean.tsv", "", set_one_name)
set_two_name <- basename(args[2])
set_two_name <- gsub("_gene-group_file.tsv", "", set_two_name)
set_two_name <- gsub("_genes_over_2sd_from_mean.tsv", "", set_two_name)

pval_tibble <- output[[2]]
as_tibble(pval_tibble) %>% 
  filter(Pval < 0.05) %>% 
  write_delim(paste0(out_path, set_one_name, set_two_name, 
                     "_significant_overlap.tsv"), 
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
