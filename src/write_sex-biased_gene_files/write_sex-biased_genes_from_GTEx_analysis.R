# this script converts ensembl gene IDs to entrez genes for 
# the genes found in GTEx analysis

library(tidyverse)

gtex <- read_delim("./signif.sbgenes.txt", delim = "\t", col_names = T)
gtex <- gtex %>% separate(gene, into = c("ensembl_gene_id", "ensembl_gene_version"), 
                          sep = "\\.")

### get ENSG to Entrez conversion ###
conversion_file_path <- "/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__ENSG-to-Entrez__All-Mappings.tsv"
gene_map <- read_delim(conversion_file_path, delim = "\t", 
                       col_names = F, col_types = cols(.default = "c"))
colnames(gene_map) <- c("ensembl_gene_id", "entrez")
gene_map <- gene_map %>% separate_rows(entrez, sep = ", ")

gtex <- left_join(gtex, gene_map, by = "ensembl_gene_id")

gtex %>% 
  filter(!is.na(entrez)) %>% 
  dplyr::select(entrez, tissue, effsize, effsize_se, lfsr) %>% 
  write_delim("./gtex_tissue_sex-biased_entrez_genes.tsv",
              delim = "\t", col_names = T)
