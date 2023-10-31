# This file will submit up to 1000 jobs at once 
library("tidyverse")

# make a directory to put qsubs and the output files
dirname <- "./sbatches_counts_agg_first_round"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (job_number in c(0:717)){
  #make sbatch
  rjobsh <- paste0("section", job_number, ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=09:00:00",
               "#SBATCH --mem=50GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=jobnum", job_number),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "",
               "",
               paste("python aggregate_counts_data_first_round.py", 
                     job_number, sep=" ")), rjobConn)
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