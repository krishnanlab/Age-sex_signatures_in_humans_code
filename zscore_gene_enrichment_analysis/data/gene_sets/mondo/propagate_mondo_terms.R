#this script filters input file, adds names from terms.tsv and writes in gmt format
#input file from mammalian phenotype ontology
tic <- as.integer(as.POSIXct(Sys.time()))

# This script will hard code what doid obo files will be used

library(tidyverse)
base_dir  <- '/mnt/research/compbio/krishnanlab/data/monarch/propagated_mondo_annotations/'

FN_direct <- paste0(base_dir,'gene_disease.9606.tsv')

#read in file as tibble with all character columns
direct_tibble <- read_delim(FN_direct, 
                            delim = "\t", 
                            col_names = TRUE,
                            col_types = cols(.default = "c")) %>% 
  select(object, subject_label)
colnames(direct_tibble) <- c("ont_id", "symbol")

#conversion file for sumbol to entrez
conversion_file <- "/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__Symbol-to-Entrez__All-Mappings.tsv"
conversion_tibble <- read_delim(conversion_file,
                                col_names = F,
                                delim = "\t",
                                col_types = cols(.default = "c"))
colnames(conversion_tibble) <- c("symbol", "gene_id")
conversion_tibble <- conversion_tibble %>% 
  separate_rows(gene_id, sep = ", ")

#convert symbol to entrez
direct_tibble <- left_join(direct_tibble, conversion_tibble, by = "symbol")
#remove rows with no entrez id, remove symbol, make unique
direct_tibble <- direct_tibble %>% 
  filter(!is.na(gene_id)) %>% 
  select(-symbol) %>% 
  distinct()

#read file with ontology terms and ancestors
FN_ancestor <- paste0(base_dir,'mondo_ancestors.tsv')
ancestor_tibble <- read_delim(FN_ancestor, delim = "\t", col_names = F)
#give colnames to ancestor tibble
colnames(ancestor_tibble) <- c("ont_id", "ancestors")
#make ancestors tibble tidy
ancestor_tibble <- separate_rows(ancestor_tibble, ancestors, sep = ", ")

# This section is the one that does the propagation
joined_tibble <- left_join(direct_tibble, ancestor_tibble, by = "ont_id")
gene_id_tibble <- joined_tibble %>% select(gene_id, ont_id)
gene_anc_tibble <- joined_tibble %>% select(gene_id, ancestors)
#rename cols of gene_anc so we can bind rows
colnames(gene_anc_tibble) <- c("gene_id", "ont_id")
#bind rows to get all gene associations - direct and ancestor relationships
joined_tibble <- bind_rows(gene_id_tibble, gene_anc_tibble)
#remove duplicated rows
joined_tibble <- distinct(joined_tibble)

#read in terms file
FN_terms <- paste0(base_dir,'mondo_terms.tsv')
terms_tibble <- read_delim(FN_terms, delim = "\t", col_names = F)
colnames(terms_tibble) <- c("ont_id", "ont_name")

# add ont names missing from that file
tmp <- read_delim(FN_direct, 
                  delim = "\t",
                  col_names = TRUE,
                  col_types = cols(.default = "c")) %>% 
  select(object, object_label) %>% 
  distinct() %>%
  filter(!object %in% terms_tibble$ont_id)
colnames(tmp) <- c("ont_id", "ont_name")
terms_tibble <- bind_rows(terms_tibble, tmp)
terms_tibble <- distinct(terms_tibble)

#add ont_name to tibble
joined_tibble <- left_join(joined_tibble, terms_tibble, by = "ont_id")

#write files
joined_tibble %>% 
  select(gene_id, ont_name) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mondo/mondo_gene-group_file.tsv",
            col_names = TRUE)

joined_tibble %>%
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mondo/mondo_gene-ont_id-ont_name_file.tsv", col_names = TRUE)
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
