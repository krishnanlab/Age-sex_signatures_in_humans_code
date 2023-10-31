# this script writes the ENST to Entrez
# conversion file using ensembl release 96
# to convert refine.bio data to Entrez 
# keeps only transcripts present in all refine.bio data
# run before aggragating data of course

tic <- as.integer(as.POSIXct(Sys.time()))
library(biomaRt)
library(tidyverse)

# get ENST to Entrez conversion
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      version = 96)
gene_map <- getBM(attributes = c("entrezgene",
                                 "ensembl_transcript_id"),
                  mart = ensembl)
gene_map <- as_tibble(gene_map)
# transcript = ENST, gene = Entrez
colnames(gene_map) <- c("gene", "transcript")

# refine.bio projects dir
data_dir <- "/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio"
# get transcripts common to all refine.bio data
# (all refine.bio datasets aligned to one of two genome verions)
# (one version = 185401 transcripts, one = 189440 transcripts)
low <- read_delim(paste0(data_dir, "/SRP048889/SRR1611777_quant.sf"), delim="\t", col_names=T) 
high <- read_delim(paste0(data_dir, "/SRP166888/SRR8112266_quant.sf"), delim="\t", col_names=T)
common_transcripts <- low %>% 
  filter(Name %in% high$Name) %>% 
  dplyr::select(Name) %>% 
  rename(transcript = Name)

gene_map <- left_join(common_transcripts, gene_map, by = "transcript")
# there are 36,423 refine.bio common_transcripts which do not map to an Entrez ID
# 149,046 do map
gene_map <- gene_map %>% 
  filter(!is.na(gene)) %>% 
  write_delim("/mnt/home/john3491/projects/age-sex-prediction/data/refine.bio/Entrez-ENST_conversion_for_refine.bio_transcripts.tsv",
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
