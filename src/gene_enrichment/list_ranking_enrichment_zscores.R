# code modified from Alex
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(parallel)
args <- commandArgs(TRUE)
#Algorithm steps
# For each gobp:
#         1. get genes in that gobp
#         2. get average score of these genes
#         3. sample "termsize" genes from all genes/scores many times, term size 
#             being size of cur gobp and get average score of each sample
#         4. Take average and sd of these averaged sample scores
#         5. Calculate z-score

#function has to be called for an arbitrarily large number of iterations.
k=100000
getavgs<-function(iter,inputgenes,termsize){
  return(mean(sample(inputgenes$scores,size=termsize)))
}

# genes with ranking to use
inputgenes = read_delim(args[1], delim = "\t", col_names = T)
colnames(inputgenes) = c("genes", "scores")

# terms to find enrichment/zscore for
# first col is gene list, second col is group or term
terms = read_tsv(args[2], col_names = T)
colnames(terms) = c("gene_id", "group")

out_path <- "~/projects/age-sex-prediction/results/age_prediction_by_sex/gene_enrichment_analysis/output_files/"
output <- tibble(group = "term", zscore = 1.2)
#for each go term you want enrichment for
for(term in unique(terms$group)){
  #get current go term and the genes for it
  curgobp=terms %>% filter(group==term)
  #curgobp=separate_rows(curgobp,gene_id,convert=T)

  #avg score of gobp genes of interest that are in original list
  gobpscores=inputgenes %>% filter(genes %in% curgobp$gene_id)
  if (dim(gobpscores)[1] == 0){
    next
  }
  avggobpgenes=mean(gobpscores$scores)

  #call function. mclapply to speed this up drastically if you have >1 core.
  sampsize=nrow(curgobp)
  iters=1:k
  avgs=mclapply(iters,getavgs,inputgenes,sampsize,mc.cores=detectCores()-1)
  avgs=unlist(avgs)
  
  #avg and sd of avgs
  avgavgs=mean(avgs)
  sdavgs=sd(avgs)
  
  #zscore
  zscore=(avggobpgenes- avgavgs)/ (sdavgs)
  term_tibble <- tibble(group = term, zscore = zscore)
  output <- bind_rows(output, term_tibble)
}

file_name <- basename(args[1])
file_name <- gsub("_model_weights.tsv", "", file_name)
term_file_name <- basename(args[2])
term_file_name <- gsub("_gene-group_file.tsv", "", term_file_name)
term_file_name <- gsub("_genes_over_2sd_from_mean.tsv", "", term_file_name)
output <- output[-1,]
output %>% 
  write_delim(paste0(out_path, file_name, "_ranks_", term_file_name, "_terms_zscores_from_list_ranking.tsv"),
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
