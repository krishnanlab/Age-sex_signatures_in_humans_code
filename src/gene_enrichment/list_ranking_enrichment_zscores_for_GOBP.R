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

#Go terms I used to test this with
gos = read_tsv("/mnt/research/compbio/krishnanlab/data/MyGeneInfo/20201029_Entrez_Multiple-Species/Gene_Ontology/propagated_annotations/GO__propagated-annotations__Homo_sapiens__Entrez__BP__EXP_IDA_IPI_IMP_IGI_TAS_IC.tsv")
gos = gos %>% filter(gene_count < 500)
#Fake gene scores for testng. I made a tibble here but all of this should easily
#translate into a named list if that is what is actually being used
allgenes = (separate_rows(gos,gene_id,convert=T))$gene_id %>% unique
inputgenes = read_delim(args[1], delim = "\t", col_names = T)
colnames(inputgenes) = c("genes", "scores")

out_path <- "~/projects/age-sex-prediction/results/age_prediction_by_sex/gene_enrichment_analysis/output_files/"
output <- tibble(gobp = "term", zscore = 1.2)
#for each go term you want enrichment for
for(term in unique(gos$go_id)){
  #get current go term and the genes for it
  curgobp=gos %>% filter(go_id==term)
  curgobp=separate_rows(curgobp,gene_id,convert=T)
  
  #avg score of gobp genes of interest that are in original list
  gobpscores=inputgenes %>% filter(genes %in% curgobp$gene_id)
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
  term_tibble <- tibble(gobp = term, zscore = zscore)
  output <- bind_rows(output, term_tibble)
}

file_name <- basename(args[1])
file_name <- gsub("_model_weights.tsv", "", file_name)
output <- output[-1,]
output %>% 
  write_delim(paste0(out_path, file_name, "_GOBP_zscores_from_list_ranking.tsv"),
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
