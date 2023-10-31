tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

path <- "/mnt/research/compbio/krishnanlab/data/GTEx/sex-biased_analyses/sex-stratified_tissue_genes/"
conversion_file <- "/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/ID_conversions/Homo_sapiens__ENSG-to-Entrez__All-Mappings.tsv"

gene_conversion <- read_tsv(conversion_file, col_types = cols(.default = "c"),
                            col_names = F)
colnames(gene_conversion) <- c("ENSG", "Entrez")
gene_conversion <- gene_conversion %>% 
  separate_rows(Entrez, sep = ",")

# read in phenotype data (SUBJID|SEX|AGE|DTHHRDY)
pheno <- read_tsv(paste0(path, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))
# read in attribute data, select sample ID, general tissue
attr <- read_tsv(paste0(path, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")) %>% 
  select(SAMPID, SMTS) %>% 
  distinct()
# gene tpm for all samples
tpm <- read_tsv(paste0(path, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"),
                skip = 2)

# sep expression by sex
female_sample_prefixes <- pheno %>% 
  filter(SEX == 2) %>% 
  pull(SUBJID) %>% 
  unique()
male_sample_prefixes <- pheno %>% 
  filter(SEX == 1) %>% 
  pull(SUBJID) %>% 
  unique()

female_tpm <- tpm %>% 
  select(Name, contains(c(female_sample_prefixes)))

male_tpm <- tpm %>% 
  select(Name, contains(c(male_sample_prefixes)))

# find median tpm expression for each tissue
# female 
out <- tibble()

# female tissues to loop through
fem_tis <- attr %>% 
  filter(str_detect(SAMPID, paste(female_sample_prefixes, collapse="|"))) %>% 
  pull(SMTS) %>% 
  unique()

for (cur_tissue in fem_tis){
  tissue_samples <- attr %>% 
    filter(SMTS == cur_tissue) %>% 
    pull(SAMPID)
    tmp <- female_tpm %>% 
      select(Name, contains(c(tissue_samples))) %>% 
      rowwise() %>% 
      mutate(median_tpm = median(c_across(where(is.numeric)))) %>% 
      select(Name, median_tpm) %>% 
      mutate(tissue = cur_tissue) %>% 
      mutate(sex = "female")
    out <- bind_rows(out, tmp)
}

# male tissues to loop through
mal_tis <- attr %>% 
  filter(str_detect(SAMPID, paste(male_sample_prefixes, collapse="|"))) %>% 
  pull(SMTS) %>% 
  unique()

# male
for (cur_tissue in mal_tis){
  tissue_samples <- attr %>% 
    filter(SMTS == cur_tissue) %>% 
    pull(SAMPID)
  tmp <- male_tpm %>% 
    select(Name, contains(c(tissue_samples))) %>% 
    rowwise() %>% 
    mutate(median_tpm = median(c_across(where(is.numeric)))) %>% 
    select(Name, median_tpm) %>% 
    mutate(tissue = cur_tissue) %>% 
    mutate(sex = "male")
  out <- bind_rows(out, tmp)
}

out <- out %>% 
  separate(Name, into = c("ENSG", "tmp"), sep = "\\.")

out <- left_join(out, gene_conversion, by = "ENSG")

out %>% 
  select(ENSG, Entrez, sex, tissue, median_tpm) %>% 
  write_tsv("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sex-strat_gtex_tissue-spec/GTEX_median_tpm_for_tissue_in_each_sex.tsv")

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))






