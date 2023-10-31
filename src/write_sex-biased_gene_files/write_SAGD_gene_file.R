# write tidy version of downloaded species files from SAGD
# with more info from Runinfo file
library(tidyverse)

# read in all files, bind into one df
species_files <- list.files("./downloaded_files", 
                            full.names = T, pattern = "\\.csv")
species_data <- lapply(species_files, read_delim, delim = ",", col_names = T)
names(species_data) <- basename(species_files)
species_data <- bind_rows(species_data, .id = "SAGD_ID")
# name unnamed column
species_data <- species_data %>% rename(ensembl_gene_id = X1)
# remove ".csv" from names
species_data$SAGD_ID <- gsub(".csv", "", species_data$SAGD_ID)

# join with run info
runinfo <- read_delim("../Runinfo.csv", delim=",", col_names=T)
runinfo <- runinfo %>% 
  dplyr::select(SAGD_ID, SRA_Project, SRA_Experiment, Tissue, Stage) %>% 
  distinct()

species_data <- left_join(species_data, runinfo, by = "SAGD_ID")

species_data <- species_data %>% 
  dplyr::select(Tissue, ensembl_gene_id, Stage,
                baseMean, log2FoldChange, lfcSE, stat, pvalue,
                padj, FPKM_M, FPKM_F, SAGD_ID, SRA_Project, SRA_Experiment) 

species_data %>% 
  write_delim("./SAGD_human_ENSG_genes_all_scores.tsv",
              delim = "\t", col_names = T)

species_data %>% 
  filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% 
  filter(padj < 0.05) %>% 
  write_delim("SAGD_human_sex-biased_ENSG_genes_L2FCover2_padj-under0.05.tsv",
              delim = "\t", col_names = T)

### get ENSG to Entrez conversion ###
conversion_file_path <- "/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__ENSG-to-Entrez__All-Mappings.tsv"
gene_map <- read_delim(conversion_file_path, delim = "\t", 
                       col_names = F, col_types = cols(.default = "c"))
colnames(gene_map) <- c("ensembl_gene_id", "entrez")
gene_map <- gene_map %>% separate_rows(entrez, sep = ", ")

species_data <- left_join(species_data, gene_map, by = "ensembl_gene_id")

species_data <- species_data %>% 
  dplyr::select(Tissue, entrez, Stage,
                baseMean, log2FoldChange, lfcSE, stat, pvalue,
                padj, FPKM_M, FPKM_F, SAGD_ID, SRA_Project, SRA_Experiment) %>% 
  filter(!is.na(log2FoldChange))

species_data %>% 
  write_delim("SAGD_human_entrez_genes_all_scores.tsv",
              delim = "\t", col_names = T)

species_data %>% 
  filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% 
  filter(padj < 0.05) %>%  
  write_delim("./SAGD_human_sex-biased_entrez_genes_L2FCover2_padj-under0.05.tsv",
              delim = "\t", col_names = T)
