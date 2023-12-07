# code modified from Alex
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(parallel)
args <- commandArgs(TRUE)
# args
# args[1] = file path for either:
#                model weights for given fold in given sex and age group
#                ranks for naive sex prediction from given age group
# args[2] = rank_method (elasticnet LR models or naive sex pred), ['LR' or 'naive']
# args[3] = data_type of file path (rnaseq or microarray)
# args[4] = gene_set, one of:
#                       disgenetC
#                       disgenetF
#                       gobp
#                       gtex
#                       guo
#                       gwas
#                       mp
#                       sagd
#                       mondo
#                       gtex_tissue-spec
#                       sex-strat_gtex_tissue-spec

# Algorithm steps
# For each gobp:
#         1. get genes in that gobp
#         2. get average score of these genes
#         3. sample "termsize" genes from all genes/scores many times, term size 
#             being size of cur gobp and get average score of each sample
#         4. Take average and sd of these averaged sample scores
#         5. Calculate z-score

#function has to be called for an arbitrarily large number of iterations.
k=100000
getavgs <- function(iter, inputgenes, termsize){
  return(mean(sample(inputgenes$scores, size = termsize)))
}

# genes with ranking to use
inputgenes = read_delim(args[1], delim = "\t", col_names = T)
colnames(inputgenes) = c("genes", "scores")

# terms to find enrichment/zscore for
# first col is gene list, second col is group or term
if (args[2] == "LR"){
  rank_method <- "elasticnet_LR_model_weights"
}
if (args[2] == "naive"){
  rank_method <- "naive_sex_prediction_ranks"
}
data_type <- args[3]
gene_set <- args[4]


if (gene_set == "disgenetC"){
  terms_path <- "../data/gene_sets/DisGeNet/disgenet_curated_gene-group_file.tsv"
}
if (gene_set == "disgenetF"){
  terms_path <- "../data/gene_sets/DisGeNet/disgenet_full_gene-group_file.tsv"
}
if (gene_set == "gobp"){
  terms_path <- "../data/gene_sets/GOBPs/GOBP_gene-group_file.tsv"
}
if (gene_set == "gtex"){
  terms_path <- "../data/gene_sets/GTEx/gtex-genes_gene-group_file.tsv"
}
if (gene_set == "guo"){
  terms_path <- "../data/gene_sets/Guo_et_al/guo_gene-group_file.tsv"
}
if (gene_set == "gwas"){
  terms_path <- "../data/gene_sets/GWAS/GWASatlas_top25_gene-group_file.tsv"
}
if (gene_set == "mp"){
  terms_path <- "../data/gene_sets/mp/mpo_mouse-human_phenotypes_gene-group_file.tsv"
}
if (gene_set == "sagd"){
  terms_path <- "../data/gene_sets/SAGD/sagd_gene-group_file.tsv"
}
if (gene_set == "mondo"){
  terms_path <- "../data/gene_sets/mondo/mondo_gene-group_file.tsv"
}
if (gene_set == "gtex_tissue-spec"){
  terms_path <- "../data/gene_sets/gtex_tissue-spec/gtex_tissue-spec_gene-group_file.tsv"
}
if (gene_set == "sctype"){
  terms_path <- "../data/gene_sets/sctype/sctype_gene-group_file.tsv"
}
if (gene_set == "sex-strat_gtex_tissue-spec"){
  terms_path <- "../data/gene_sets/sex-strat_gtex_tissue-spec/sex-strat_gtex_tissue-spec_gene-group_file.tsv"
}

terms = read_tsv(terms_path, col_names = T)
colnames(terms) = c("gene_id", "group") 

out_path <- paste0("../results/",
                   rank_method, "/", data_type, "/", gene_set, "/")

output <- tibble(group = "term", zscore = 1.22)
#for each go term you want enrichment for
for(term in unique(terms$group)){
  #get current go term and the genes for it
  curterm=terms %>% filter(group==term)
  
  #avg score of gobp genes of interest that are in original list
  termscores=inputgenes %>% filter(genes %in% curterm$gene_id)
  if (nrow(termscores) < 10){
    next
  }
  avgtermgenes=mean(termscores$scores)
  
  #call function. mclapply to speed this up drastically if you have >1 core.
  sampsize=nrow(termscores)
  iters=1:k
  avgs=mclapply(iters,getavgs,inputgenes,sampsize,mc.cores=detectCores()-1)
  avgs=unlist(avgs)
  
  #avg and sd of avgs
  avgavgs=mean(avgs)
  sdavgs=sd(avgs)
  
  #zscore
  zscore=(avgtermgenes- avgavgs)/ (sdavgs)
  term_tibble <- tibble(group = term, zscore = zscore)
  output <- bind_rows(output, term_tibble)
}

file_name <- basename(args[1])
file_name <- gsub("_model_weights.tsv", "", file_name)
file_name <- gsub("_signed_gene_cut_data_all_chromosomes.tsv", "", file_name)

output <- output[-1,]
output %>% 
  write_delim(paste0(out_path, file_name, "___ranks_", gene_set, "_terms_zscores_from_list_ranking.tsv"),
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
