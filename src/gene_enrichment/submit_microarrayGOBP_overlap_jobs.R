tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

file_list <- list.files("~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses", 
                        pattern = "*_model_weights.tsv", full.names = T)

# make a directory to put sbatches
dirname <- paste0("./sbatches_GOBP_microarray_overlap")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

args = c("1 10", "11 20", "21 30", "31 40", "41 50",
         "51 60", "61 70", "71 80", "81 90", "91 100",
         "101 110", "111 120", "121 130", "131 140", "141 150",
         "151 160", "161 170", "171 180", "181 190", "191 200",
         "200 210", "210 216")

for (arg in args) {
  #make sbatch
  rjobsh <- paste0("GOBP_overlap", gsub(" ","_",arg), ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=4:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=microarray_GOBPoverlap_", gsub(" ","_",arg)),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "module purge",
               "module load GCC/8.3.0 OpenMPI/3.1.4",
               "module load R/4.0.2",
               "",
               paste("Rscript GOBP_enrichment_for_microarray_models.R",
                     arg, sep = " "),
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


