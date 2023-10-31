#this script filters input file, adds names from terms.tsv and writes in gmt format
#input file from mammalian phenotype ontology
tic <- as.integer(as.POSIXct(Sys.time()))

# This script will hard code what doid obo files will be used

library(tidyverse)
base_dir  <- '/mnt/research/compbio/krishnanlab/data/ontologies/mammalian_phenotype_ontology/'

FN_direct <- paste0(base_dir,'MGI_GenePheno.rpt')

#read in file as tibble with all character columns
#it reads in a 6th column which is all NA, I think induced by trailing tabs
#but we only need the human entrez column (2) and the MP annotations (5)
#get rid of genes/rows that have no annotations
direct_tibble <- read_delim(FN_direct, 
                            delim = "\t", 
                            col_names = FALSE,
                            col_types = cols(.default = "c")) %>% 
  select(X5, X7)
colnames(direct_tibble) <- c("ont_id", "gene_id")

#read in terms file
FN_terms <- paste0(base_dir,'mp_terms.tsv')
terms_tibble <- read_delim(FN_terms, delim = "\t", col_names = F)
colnames(terms_tibble) <- c("ont_id", "ont_name")

#read file with ontology terms and ancestors
FN_ancestor <- paste0(base_dir,'mp_ancestors.tsv')
ancestor_tibble <- read_delim(FN_ancestor, delim = "\t", col_names = F)
#give colnames to ancestor tibble
colnames(ancestor_tibble) <- c("ont_id", "ancestors")
#make ancestors tibble tidy
ancestor_tibble <- separate_rows(ancestor_tibble, ancestors, sep = ", ")
direct_tibble <- direct_tibble[!duplicated(direct_tibble),]
direct_tibble <- separate_rows(direct_tibble, ont_id, sep = ", ")

# This section is the one that does the propagation
joined_tibble <- left_join(direct_tibble, ancestor_tibble, by = "ont_id")
gene_id_tibble <- joined_tibble %>% select(gene_id, ont_id)
gene_anc_tibble <- joined_tibble %>% select(gene_id, ancestors)
#rename cols of gene_anc so we can bind rows
colnames(gene_anc_tibble) <- c("gene_id", "ont_id")
#bind rows to get all gene associations - direct and ancestor relationships
joined_tibble <- bind_rows(gene_id_tibble, gene_anc_tibble)
#remove duplicated rows
joined_tibble <- joined_tibble[!duplicated(joined_tibble),]

#add ont_name to tibble
joined_tibble <- left_join(joined_tibble, terms_tibble, by = "ont_id")

#convert mouse genes to human
hom_file <- read_delim(paste0(base_dir, "HOM_MouseHumanSequence.rpt"), 
                              col_names = T, delim = "\t") %>% 
  select(`DB Class Key`, `NCBI Taxon ID`, `EntrezGene ID`, `Mouse MGI ID`)
#make mouse and human tables to join by `DB Class Key`
mouse <- hom_file %>% 
  filter(`NCBI Taxon ID` == 10090) %>% 
  select(`DB Class Key`, `Mouse MGI ID`)
human <- hom_file %>% 
  filter(`NCBI Taxon ID` == 9606) %>% 
  select(`DB Class Key`, `EntrezGene ID`)
#join, get conversion table
gene_conversion <- inner_join(human, mouse, by = "DB Class Key")
# mouse gene col named "gene_id" so it can be joined with other tibble to convert
colnames(gene_conversion) <- c("key", "entrez", "gene_id")
#convert
joined_tibble <- left_join(joined_tibble, gene_conversion, by = "gene_id")
#select desired cols
joined_tibble <- joined_tibble %>% 
  select(entrez, ont_id, ont_name, gene_id)
#change col names to be more informative
colnames(joined_tibble) <-  c("human_entrez", "mpo_id", "mpo_name", "mouse_gene")

#write files
joined_tibble %>% 
  select(human_entrez, mpo_name) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mp/mpo_mouse-human_phenotypes_gene-group_file.tsv",
            col_names = TRUE)

joined_tibble %>%
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mp/mpo_mouse-human_gene-mpoID-mpoName-entrez_file.tsv", col_names = TRUE)
#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
