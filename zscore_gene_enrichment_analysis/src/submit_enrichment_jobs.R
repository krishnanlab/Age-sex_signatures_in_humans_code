tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
args <- commandArgs(TRUE)
# args:
# args[1] = rank_method (elasticnet LR models or naive sex pred), ['LR' or 'naive']
# args[2] = data_type (data type used in ranking method), ['rnaseq' or 'microarray']
# args[3] = gene_set = one of:
#                       'disgenetC'
#                       'disgenetF'
#                       'gobp'
#                       'gtex'
#                       'guo'
#                       'gwas'
#                       'mp'
#                       'sagd'
#                       'mondo'
#                       'gtex_tissue-spec'
#                       'sex-strat_gtex_tissue-spec'

rank_method <- args[1]
data_type <- args[2]
gene_set <- args[3]

if (rank_method == "LR"){
  file_list <- list.files(paste0("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/elasticnet_LR_model_weights/",
                                 data_type), 
                          pattern = "*_model_weights.tsv", full.names = T)
}

if (rank_method == "naive"){
  file_list <- list.files(paste0("~/projects/age-sex-prediction/final_list-rank-zscore_gene_enrichment_analyses/data/naive_sex_prediction_signed_ranks/",
                                 data_type), 
                          pattern = "*_signed_gene_cut_data_all_chromosomes.tsv", full.names = T)
}


# make a directory to put sbatches
dirname <- paste0("./sbatches_gene_list_ranking_zscores_", rank_method, "_", data_type, "_", gene_set)
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (file_name in file_list) {
  #make sbatch
  rjobsh <- paste0("list_ranking", basename(file_name), ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=4:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=25",
               paste0("#SBATCH --job-name=term_overlap_", basename(file_name)),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "module purge",
               "module load GCC/8.3.0 OpenMPI/3.1.4",
               "module load R/4.0.2",
               "",
               paste("Rscript list_ranking_enrichment_zscores.R",
                     file_name, rank_method, data_type, gene_set, sep = " "),
               ""), rjobConn)
  close(rjobConn)
  
  system(paste0("sbatch ", paste0(dirname,"/",rjobsh)))
  
  njobs <- system("squeue -u  john3491 | wc -l", intern=TRUE)
  njobs <- as.numeric(njobs)
  while(njobs > 1000) {
    Sys.sleep(360)
    njobs <- system("squeue -u  john3491 | wc -l", intern=TRUE)
    njobs <- as.numeric(njobs)
  }
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

