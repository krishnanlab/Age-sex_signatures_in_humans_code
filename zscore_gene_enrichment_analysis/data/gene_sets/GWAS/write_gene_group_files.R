tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(qvalue)

#read in pvalues 
p_vals <- read_tsv("/mnt/research/compbio/krishnanlab/tmp/to_Kayla_from_Alex/magma.P.r3.txt")
meta <- read_tsv("/mnt/research/compbio/krishnanlab/tmp/to_Kayla_from_Alex/gwasATLAS_v20191115.txt") %>% 
  select(id, Trait)

gene_conversion <- read_tsv("/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__ENSG-to-Entrez__All-Mappings.tsv",
                            col_names = F, col_types = cols(.default = "c")) %>% 
  separate_rows(X2, sep = ", ")
colnames(gene_conversion) <- c("ensg", "entrez")

# outputs for FDR = 0.05, 0.1
output5 <- tibble(id = 1, genes = "ENSG22")
output10 <- tibble(id = 1, genes = "ENSG22")

for (n in 2:4757){
  #select one trait at a time
  p <- p_vals %>% 
    select(c(1, all_of(n)))
  colnames(p) <- c("ensg", "pvalues")
  p <- p %>% filter(!is.na(pvalues))
  p_vec <- pull(p, pvalues)
  # try to calculate q values
  # if not possible (pi0 cannot be estimated from the data),
  # force BH procedure to be used by setting pi0 = 1
  q_vec <- tryCatch({qvalue::qvalue(p_vec)$qvalues},
           error = function(e) {return(qvalue::qvalue(p_vec, pi0 = 1)$qvalues)})
  # make ensg/pvalues/qvalues tbl
  tbl <- bind_cols(p, q_vec)
  colnames(tbl) <- c("ensg", "pvalues", "qvalues")
  # pull genes under FDR 0.1
  tbl <- tbl %>% 
    filter(qvalues < 0.1)
  tbl <- left_join(tbl, gene_conversion, by = "ensg")
  trait_genes <- tbl %>% 
    filter(!is.na(entrez)) %>% 
    pull(entrez)
  trait_out <- tibble(id = n, genes = trait_genes)
  output10 <- bind_rows(output10, trait_out)
  trait_genes <- tbl %>% 
    filter(qvalues < 0.05) %>% 
    pull(entrez)
  trait_out <- tibble(id = n, genes = trait_genes)
  output5 <- bind_rows(output5, trait_out)
}
# get rid of initializing line
output5 <- output5[-1,]
output10 <- output10[-1,]
# adjust id so it reflects trait number
# (because first column = genes, trait 1 is col 2 and so on)
output5 <- output5 %>% mutate(id = id - 1)
output10 <- output10 %>% mutate(id = id - 1)

# FDR = 0.05 (output5) number genes:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.0     2.0     7.0   110.1    54.0 10235.0

# FDR = 0.1 (output10) number genes:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.0     2.0     9.0   149.1    60.0 14602.0

# output for top 25, 50 genes
output25 <- tibble(id = 1, genes = "ENSG22")
output50 <- tibble(id = 1, genes = "ENSG22")

for (n in 2:4757){
  p <- p_vals %>% 
    select(c(1, all_of(n)))
  colnames(p) <- c("ensg", "pvalues")
  p <- left_join(p, gene_conversion, by = "ensg")
  # NAs automatically sort to end w/ arrange
  p <- p %>% 
    filter(!is.na(entrez)) %>% 
    arrange(pvalues)
  p <- p[1:50,]
  trait_genes <- p %>% 
    pull(entrez)
  trait_out <- tibble(id = n, genes = trait_genes)
  output50 <- bind_rows(output50, trait_out)
  p <- p[1:25,]
  trait_genes <- p %>% 
    pull(entrez)
  trait_out <- tibble(id = n, genes = trait_genes)
  output25 <- bind_rows(output25, trait_out)
}
# get rid of initializing line
output25 <- output25[-1,]
output50 <- output50[-1,]
# adjust n so it reflects trait number
# (because first column = genes, trait 1 is col 2 and so on)
output25 <- output25 %>% mutate(id = id - 1)
output50 <- output50 %>% mutate(id = id - 1)

# write files with trait names
output5 <- left_join(output5, meta, by = "id")
output10 <- left_join(output10, meta, by = "id")
output25 <- left_join(output25, meta, by = "id")
output50 <- left_join(output50, meta, by = "id")

# paste is and trait name
# select genes and that column
# write file
output5 %>% 
  mutate(trait = paste(id, Trait, sep = "_")) %>% 
  select(genes, trait) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GWAS/GWASatlas_FDR0.05_gene-group_file.tsv")

output10 %>% 
  mutate(trait = paste(id, Trait, sep = "_")) %>% 
  select(genes, trait) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GWAS/GWASatlas_FDR0.1_gene-group_file.tsv")

output25 %>% 
  mutate(trait = paste(id, Trait, sep = "_")) %>% 
  select(genes, trait) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GWAS/GWASatlas_top25_gene-group_file.tsv")

output50 %>% 
  mutate(trait = paste(id, Trait, sep = "_")) %>% 
  select(genes, trait) %>% 
  distinct() %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GWAS/GWASatlas_top50_gene-group_file.tsv")

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
