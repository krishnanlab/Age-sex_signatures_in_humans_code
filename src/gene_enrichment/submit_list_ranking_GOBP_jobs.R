tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

file_list <- list.files("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses", 
                        pattern = "*_model_weights.tsv", full.names = T)

# make a directory to put sbatches
dirname <- paste0("./sbatches_GOBP_list_ranking_zscores")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (file_name in file_list) {
  #make sbatch
  rjobsh <- paste0("GOBP_list_ranking", basename(file_name), ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=4:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=25",
               paste0("#SBATCH --job-name=rnaseq_overlap_", basename(file_name)),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "module purge",
               "module load GCC/8.3.0 OpenMPI/3.1.4",
               "module load R/4.0.2",
               "",
               paste("Rscript list_ranking_enrichment_zscores_for_GOBP.R",
                     file_name, sep = " "),
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

