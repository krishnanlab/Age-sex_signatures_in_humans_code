library(tidyverse)

# read in data
gobp <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GOBPs/GOBP_gene-group_file.tsv")
disc <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_curated_gene-group_file.tsv")
disf <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/DisGeNet/disgenet_full_gene-group_file.tsv")
gtex <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GTEx/gtex-genes_gene-group_file.tsv")
sagd <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/SAGD/sagd_gene-group_file.tsv")
guo <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/Guo_et_al/guo_gene-group_file.tsv")
gwas <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/GWAS/GWASatlas_top25_gene-group_file.tsv")
mp <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mp/mpo_mouse-human_phenotypes_gene-group_file.tsv")
mondo <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/mondo/mondo_gene-group_file.tsv")
ts <- read_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/gtex_tissue-spec/gtex_tissue-spec_gene-group_file.tsv")

common_genes <- read_tsv("~/projects/age-sex-prediction/data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt",
                         col_names = F) %>% pull(1)

# filter to only genes in RNAseq+microarray data
gobp <- gobp %>% filter(gene_id %in% common_genes)
disc <- disc %>% filter(geneId %in% common_genes)
disf <- disf %>% filter(geneId %in% common_genes)
gtex <- gtex %>% filter(entrez %in% common_genes)
sagd <- sagd %>% filter(entrez %in% common_genes)
guo <- guo %>% filter(entrez %in% common_genes)
gwas <- gwas %>% filter(genes %in% common_genes)
mp <- mp %>% filter(human_entrez %in% common_genes)
mondo <- mondo %>% filter(gene_id %in% common_genes)
ts <- ts %>% filter(Entrez %in% common_genes)

# filter out terms that don't have at least 10 genes
gobp_terms <- gobp %>% group_by(go_id) %>% tally() %>% filter(n > 9) %>% pull(go_id) %>% unique()
disc_terms <- disc %>% group_by(diseaseName) %>% tally() %>% filter(n > 9) %>% pull(diseaseName) %>% unique()
disf_terms <- disf %>% group_by(diseaseName) %>% tally() %>% filter(n > 9) %>% pull(diseaseName) %>% unique()
gtex_terms <- gtex %>% group_by(group) %>% tally() %>% filter(n > 9) %>% pull(group) %>% unique()
sagd_terms <- sagd %>% group_by(group) %>% tally() %>% filter(n > 9) %>% pull(group) %>% unique()
guo_terms <- guo %>% group_by(group) %>% tally() %>% filter(n > 9) %>% pull(group) %>% unique()
gwas_terms <- gwas %>% group_by(trait) %>% tally() %>% filter(n > 9) %>% pull(trait) %>% unique()
mp_terms <- mp %>% group_by(mpo_name) %>% tally() %>% filter(n > 9) %>% pull(mpo_name) %>% unique()
mondo_terms <- mondo %>% group_by(ont_name) %>% tally() %>% filter(n > 9) %>% pull(ont_name) %>% unique()
ts_terms <- ts %>% group_by(tissue) %>% tally() %>% filter(n > 9) %>% pull(tissue) %>% unique()

# filter to only terms that were used (have at least ten genes in the common set)
gobp <- gobp %>% filter(go_id %in% gobp_terms)
disc <- disc %>% filter(diseaseName %in% disc_terms)
disf <- disf %>% filter(diseaseName %in% disf_terms)
gtex <- gtex %>% filter(group %in% gtex_terms)
sagd <- sagd %>% filter(group %in% sagd_terms)
guo <- guo %>% filter(group %in% guo_terms)
gwas <- gwas %>% filter(trait %in% gwas_terms)
mp <- mp %>% filter(mpo_name %in% mp_terms)
mondo <- mondo %>% filter(ont_name %in% mondo_terms)
ts <- ts %>% filter(tissue %in% ts_terms)

output <- tibble(geneset = c("GOBP", "DisGeNetC", "DisGeNetF", "GTEx", "SAGD", "Guo", "GWAS", "MPO", "MONDO", "GTEx_tissue-spec"),
                 total_genes = c(length(unique(gobp$gene_id)),
                                 length(unique(disc$geneId)),
                                 length(unique(disf$geneId)),
                                 length(unique(gtex$entrez)),
                                 length(unique(sagd$entrez)),
                                 length(unique(guo$entrez)),
                                 length(unique(gwas$genes)),
                                 length(unique(mp$human_entrez)),
                                 length(unique(mondo$gene_id)),
                                 length(unique(ts$Entrez))),
                 number_terms = c(length(gobp_terms),
                                  length(disc_terms),
                                  length(disf_terms),
                                  length(gtex_terms),
                                  length(sagd_terms),
                                  length(guo_terms),
                                  length(gwas_terms),
                                  length(mp_terms),
                                  length(mondo_terms),
                                  length(ts_terms)))

output %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/number_genes_in_RNAseq-Mircroarray_commonset.tsv")


