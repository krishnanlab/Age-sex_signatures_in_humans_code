# convert symbol to entrez from guo analysis
library(tidyverse)

guo <- read_delim("./guo_supp.tsv", delim = "\t", col_names = T)

### get symbol to Entrez conversion ###
conversion_file_path <- "/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__Symbol-to-Entrez__All-Mappings.tsv"
gene_map <- read_delim(conversion_file_path, delim = "\t", 
                       col_names = F, col_types = cols(.default = "c"))
colnames(gene_map) <- c("hgnc_symbol", "entrez")
gene_map <- gene_map %>% separate_rows(entrez, sep = ", ")

guo <- left_join(guo, gene_map, by = "hgnc_symbol")

guo %>% 
  filter(!is.na(entrez)) %>% 
  dplyr::select(entrez, chromosome, sex_bias, tissue) %>% 
  write_delim("./guo_et_al_entrez_sex-baised_genes.tsv",
              delim = "\t", col_names = T)
