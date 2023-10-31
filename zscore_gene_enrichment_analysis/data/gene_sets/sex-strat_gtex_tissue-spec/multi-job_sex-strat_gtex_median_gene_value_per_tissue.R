tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)

sx <- args[1]
cur_tissue <- args[2]

# fix two word tissues
if (cur_tissue == "Adipose"){
  cur_tissue <- "Adipose Tissue"
}
if (cur_tissue == "Salivary"){
  cur_tissue <- "Salivary Gland"
}
if (cur_tissue == "Adrenal"){
  cur_tissue <- "Adrenal Gland"
}
if (cur_tissue == "Small"){
  cur_tissue <- "Small Intestine"
}
if (cur_tissue == "Cervix"){
  cur_tissue <- "Cervix Uteri"
}
if (cur_tissue == "Fallopian"){
  cur_tissue <- "Fallopian Tube"
}
if (cur_tissue == "Bone"){
  cur_tissue <- "Bone Marrow"
}
if (cur_tissue == "BV"){
  cur_tissue <- "Blood Vessel"
}


print(paste("args:", sx, cur_tissue, sep = " "))

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
if (sx == "female"){
  sample_prefixes <- pheno %>% 
    filter(SEX == 2) %>% 
    pull(SUBJID) %>% 
    unique()
} else {
  sample_prefixes <- pheno %>% 
    filter(SEX == 1) %>% 
    pull(SUBJID) %>% 
    unique()
}

tpm <- tpm %>% 
  select(Name, contains(c(sample_prefixes)))

# find median tpm expression for tissue in correct sex

tissue_samples <- attr %>% 
    filter(SMTS == cur_tissue) %>% 
    pull(SAMPID)
out <- tpm %>% 
    select(Name, contains(c(tissue_samples))) %>% 
    rowwise() %>% 
    mutate(median_tpm = median(c_across(where(is.numeric)))) %>% 
    select(Name, median_tpm) %>% 
    mutate(tissue = cur_tissue) %>% 
    mutate(sex = sx)

out <- out %>% 
  separate(Name, into = c("ENSG", "tmp"), sep = "\\.")

out <- left_join(out, gene_conversion, by = "ENSG")

out %>% 
  select(ENSG, Entrez, sex, tissue, median_tpm) %>% 
  write_tsv(paste0("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/gene_sets/sex-strat_gtex_tissue-spec/ind_sex-tissue_outputs/GTEX_median_tpm_for_", sx, "_", gsub(" ", "_", cur_tissue), ".tsv"))

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))




