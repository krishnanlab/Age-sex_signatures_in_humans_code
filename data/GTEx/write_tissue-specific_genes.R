tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

# Function requires a vector with expression of one gene in different tissues.
fTau <- function(x){
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      # print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    # print("No data for this gene available.")
  } 
  return(res)
}

# Function requires a vector with expression of one gene in different tissues.
robust_z <- function(x){
  med <- median(x)
  rz <- (x - med)/(median(abs(x-med)) * 1.4826)
  return(rz)
}

gtex <- read_delim("~/projects/age-sex-prediction/data/GTEx/Entrez_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv",
                   delim = "\t", col_names = T, col_types = cols(.default = "d", Entrez = "i"))


general_tissue_map <- tibble(spec_tissues = c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", 
                                              "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", 
                                              "Bladder", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", 
                                              "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
                                              "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", 
                                              "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", 
                                              "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", 
                                              "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", 
                                              "Cells_EBV-transformed_lymphocytes", "Cervix_Ectocervix", "Cervix_Endocervix", 
                                              "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", 
                                              "Esophagus_Mucosa", "Esophagus_Muscularis", "Fallopian_Tube", 
                                              "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", 
                                              "Kidney_Medulla", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", 
                                              "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", 
                                              "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
                                              "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", 
                                              "Thyroid", "Uterus", "Vagina", "Whole_Blood"),
                             general_tissues = c("adipose", "adipose",
                                                 "adrenal_gland", "artery", "artery", "artery",
                                                 "bladder", "brain", "brain",
                                                 "brain", "brain", 
                                                 "brain", "brain", "brain",
                                                 "brain",  "brain",  "brain",
                                                 "brain", "brain", 
                                                 "brain",  "breast", "fibroblast_cell_line",
                                                 "lymphocyte_cell_line", "cervix", "cervix",
                                                 "colon", "colon", "esophagus",
                                                 "esophagus", "esophagus", "fallopian_tube",
                                                 "heart", "heart", "kidney",
                                                 "kidney", "liver", "lung", "salivary_gland", "muscle",
                                                 "nerve", "ovary", "pancreas", "pituitary", "prostate",
                                                 "skin", "skin",
                                                 "small_intestine", "spleen", "stomach", "testis",
                                                 "thyroid", "uterus", "vagina", "blood"))
genes <- gtex$Entrez
ggtex <- gtex %>%
  select(-Entrez) %>%
  t() %>%
  as_tibble()
ggtex$tissue <- general_tissue_map$general_tissues

ggtex <- ggtex %>%
  group_by(tissue) %>%
  summarise_all(mean)

ggtex_tissues <- ggtex$tissue
ggtex <- ggtex %>%
  select(-tissue) %>%
  t() %>%
  as_tibble()
colnames(ggtex) <- ggtex_tissues
ggtex$Entrez <- genes
ggtex <- ggtex %>%
  select(Entrez, everything())

# ggtex %>%
#   write_delim("~/projects/age-sex-prediction/data/GTEx/general_tissue_Entrez_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv",
#               delim = "\t", col_names = T)

ggenes <- ggtex$Entrez
ggtex <- asinh(ggtex[,2:33])
ggtex$Entrez <- ggenes
ggtex$Entrez <- as.character(ggtex$Entrez)
# 24,814 to 20,598 genes
ggtex <- ggtex %>% 
  filter_if(is.numeric, any_vars(. >= 1))

ggenes <- ggtex$Entrez
ggtex <- ggtex %>% select(-Entrez)

gts_scores <- apply(ggtex, 1, fTau)

zggtex <- apply(ggtex, 1, robust_z)
zggtex <- t(zggtex) %>% as_tibble()
zggtex$Entrez <- ggenes
zggtex$tau_score <- gts_scores

tidy_zggtex <- gather(zggtex, key = tissue, value = robust_zscore, c(1:32))
#zgts_scores <- apply(zggtex, 1, fTau)
#tidy_ggtex <- gather(ggtex, key = tissue, value = transformed_tpm)

tidy_zggtex %>% 
  write_delim("./GTEx_expression_Entrez-tau_score-tissue-gene_zscore.tsv",
              delim = "\t", col_names = T)


tidy_zggtex <- tidy_zggtex %>% 
  filter(tau_score >= 0.8)

for (tis in unique(tidy_zggtex$tissue)){
  tis_genes <- tidy_zggtex %>% 
    filter(tissue == tis) %>% 
    filter(robust_zscore >= 2) %>% 
    pull(Entrez)
  
  tis_genes %>% 
    as.data.frame() %>% 
    write_delim(paste0("./", tis, "-specific_genes.txt"),
                delim = "\t", col_names = F)
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

