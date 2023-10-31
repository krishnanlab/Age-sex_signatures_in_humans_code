tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

age_groups <- c("infant","young_child","child","adolescent","young_adult","adult","middle_adult","older_adult","old_adult")

# make a directory to put sbatches
dirname <- paste0("./sbatches_blood_age_prediction_by_sex_elasticnet")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (sex in c("female", "male")){
  for (dt in c("rnaseq", "microarray")){
    for (ag in age_groups) {
      #make sbatch
      rjobsh <- paste0("el_", sex, "_", ag, "_", dt, "_no_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
        
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                    "#SBATCH --time=12:00:00",
                    "#SBATCH --mem=60GB",
                    "#SBATCH --nodes=1",
                    "#SBATCH --cpus-per-task=1",
                    paste0("#SBATCH --job-name=", sex, "_", ag, "_", dt, "_no_stdsclng"),
                    "#SBATCH --account=wang-krishnan",
                    "",
                    "cd $SLURM_SUBMIT_DIR",
                    "source /mnt/home/john3491/.bashrc",
                    "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                    "which python",
                    "",
                    paste0("python age_prediction_in_blood_by_sex_for_one_age_group.py --sex=",
                           sex, " --pos_age_group=", ag, " --data_type=", dt),
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
  }
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
