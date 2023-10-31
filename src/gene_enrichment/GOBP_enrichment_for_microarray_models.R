library(tidyverse)
args <- commandArgs(TRUE)
tic <- as.integer(as.POSIXct(Sys.time()))

out_path <- "~/projects/age-sex-prediction/results/age_prediction_by_sex/gene_enrichment_analysis/output_files/"

source("findEnrichedGOBP.R")

common_genes = read_delim("~/projects/age-sex-prediction/data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", 
                          delim = "\t", col_names = F) %>% 
  pull(1)
microarray_genes <- read_delim("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/microarray_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv",
                               delim = "\t", col_names = T, col_types = cols(.default = "c"))
colnames(microarray_genes) <- c("Gene", "Group")

output <- tibble(GO.ID = "GO:0045630", Term = "positive regulation",
                 Annotated = 3, Significant = 3, Expected = 0.23, 
                 pval = "0.05", FDR = 0.05, OverlappingGenes = "114548, 5590, 942", 
                 model_group = "pos_female")

all_groups <-  unique(microarray_genes$Group)
all_groups <- all_groups[as.numeric(args[1]):as.numeric(args[2])]

for (group in unique(microarray_genes$Group)){
  genes_of_interest <- microarray_genes %>% 
    filter(Group == group) %>% 
    pull(Gene)
  group_df <- findEnrichedGOBP(goi = genes_of_interest,
                               background = common_genes,
                               org.db = "org.Hs.eg.db",
                               id_type = "entrez",
                               min_size = 5,
                               max_size = 200)
  group_df <- group_df %>% 
    filter(as.numeric(pval) < 0.05) %>% 
    mutate(model_group = group) %>% 
    as_tibble()
  output <- bind_rows(output, group_df)
}
output <- output[-1,]

output %>% 
  write_delim(paste0(out_path, as.character(args[1]), as.character(args[2]), "microarray_elasticnet_age_group_models_significant_GOBP_enrichment.tsv"),
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
