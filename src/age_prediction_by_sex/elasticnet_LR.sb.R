tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

age_groups <- c("fetus","infant","young_child","child","adolescent","young_adult","adult","middle_adult","older_adult","old_adult","elderly")

src_dir <- "/mnt/home/john3491/projects/age-sex-prediction/src/age_prediction_by_sex/"

# make a directory to put sbatches
dirname <- paste0("./sbatches_age_prediction_by_sex_elasticnet")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (sex in c("female", "male")){
  for (dt in c("rnaseq", "microarray")){
    for (ss in c("with", "without")){
      for (ag in age_groups) {
    #make sbatch
    rjobsh <- paste0("el_", sex, "_", ag, "_", dt, "_", ss, "_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
    
    rjobConn <- file(paste0(dirname,"/",rjobsh))
    writeLines(c("#!/bin/sh -login",
                 "#SBATCH --time=12:00:00",
                 "#SBATCH --mem=60GB",
                 "#SBATCH --nodes=1",
                 "#SBATCH --cpus-per-task=1",
                 paste0("#SBATCH --job-name=", sex, "_", ag, "_", dt, "_", ss, "_stdsclng"),
                 "#SBATCH --account=wang-krishnan",
                 "",
                 "cd $SLURM_SUBMIT_DIR",
                 "source /mnt/home/john3491/.bashrc",
                 "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                 "which python",
                 "",
                 paste0("python LR_elasticnet_3F_CV_for_one_age_group.py --sex=",
                        sex, " --pos_age_group=", ag,
                        " --data_type=", dt, " --std_scale=", ss),
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
}

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
