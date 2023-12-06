library(tidyverse)

rf <- read_tsv("~/projects/age-sex-prediction/blood-only_analysis/data/naive_sex_prediction_signed_ranks/rnaseq_fine_age_groups_final_gene_cut_data_fixed_step_all_chromosomes.tsv")
mf <- read_tsv("~/projects/age-sex-prediction/blood-only_analysis/data/naive_sex_prediction_signed_ranks/microarray_fine_age_groups_final_gene_cut_data_fixed_step_all_chromosomes.tsv")

common_genes <- read_tsv("~/projects/age-sex-prediction/data/refine.bio/common_Entrez_IDs_gpl570-refine.bio.txt", col_names = F) %>% 
  pull(1)

fine_ags <- c("infant", "young_child", "child", "adolescent", "young_adult", 
              "adult", "middle_adult", "older_adult", "old_adult")

rf <- rf %>% 
  filter(gene %in% common_genes)

mf <- mf %>% 
  filter(gene %in% common_genes)

no_var_genes <- setdiff(common_genes, unique(rf$gene))

rfaddition <-  tibble(age_group = rep(fine_ags, each = 164), gene = rep(no_var_genes, times = 9), cutoff = 'none', 
                      balanced_accuracy = 0.5, higher_sex = 'none')

rf <- bind_rows(rf, rfaddition)

################# add gene symbols from other files #########################
gene_sym <- read_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/20step_naive_sex_prediction/rnaseq_signed_balanced_accuracy_data.tsv",
                     col_names = T) %>% 
  select(gene, symbol) %>% 
  distinct()

# choose one symbol for duplicate genes with diff symbols
gene_sym <- gene_sym %>% 
  filter(!(gene == 1394 & symbol == "LINC02210-CRHR1")) %>% 
  filter(!(gene == 23596 & symbol == "CHML")) %>% 
  filter(!(gene == 3117 & symbol == "HLA-DQA2")) %>% 
  filter(!(gene == 3802 & symbol == "KIR2DS1")) %>% 
  filter(!(gene == 414243 & symbol == "LINC00856")) %>% 
  filter(!(gene == 100132596 & symbol == "XG"))

# remove genes that have both an NA line and a line with a symbol
tmp <- gene_sym$gene
tmp <- tmp[duplicated(tmp)]
gene_sym <- gene_sym %>% 
  filter(!(gene %in% tmp & is.na(symbol)))

# add symbols
rf <- left_join(rf, gene_sym, by = "gene")
mf <- left_join(mf, gene_sym, by = "gene")

### add chromosomes
gene_chr <- read_tsv("~/projects/age-sex-prediction/data/gene_chromosomes/clean_mygeneinfo_homo_sapiens_entrez_genes_chromosome_locations.tsv") %>% 
  select(-map_location)

rf <- left_join(rf, gene_chr, by = "gene")
mf <- left_join(mf, gene_chr, by = "gene")

rf %>% 
  write_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/blood_small_step_naive_sex_prediction/blood_rnaseq_unsigned_balanced_accuracy_info_file.tsv")

mf %>% 
  write_tsv("~/projects/age-sex-prediction/paper_figures/final_results_files/blood_small_step_naive_sex_prediction/blood_microarray_unsigned_balanced_accuracy_info_file.tsv")



