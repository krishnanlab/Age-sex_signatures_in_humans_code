tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
# this submits all of the no-asinh transformation rnaseq jobs  for l2/elasticnet and fine/coarse age groups

fine_age_groups <- c("fetus","infant","young_child","child","adolescent","young_adult","adult","middle_adult","older_adult","old_adult","elderly")
coarse_age_groups <- c("fetus","infant","child","adolescent","adult","older_adult","elderly")

src_dir <- "/mnt/home/john3491/projects/age-sex-prediction/src/age_prediction_by_sex/"

# make a directory to put sbatches
dirname <- paste0("./sbatches_age_prediction_by_sex_rnaseq_no-asinh")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

# elastic net fine age groups 
for (sex in c("female", "male")){
  for (ss in c("with", "without")){
    for (ag in fine_age_groups) {
      #make sbatch
      rjobsh <- paste0("en_", sex, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
      
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                   "#SBATCH --time=12:00:00",
                   "#SBATCH --mem=60GB",
                   "#SBATCH --nodes=1",
                   "#SBATCH --cpus-per-task=1",
                   paste0("#SBATCH --job-name=", "en_", sex, "_", ag, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng"),
                   "#SBATCH --account=wang-krishnan",
                   "",
                   "cd $SLURM_SUBMIT_DIR",
                   "source /mnt/home/john3491/.bashrc",
                   "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                   "which python",
                   "",
                   paste0("python LR_elasticnet_3F_CV_for_one_age_group.py --sex=",
                          sex, " --age_group_type=fine", " --pos_age_group=", ag,
                          " --data_type=rnaseq", " --std_scale=", ss, " --no-asinh"),
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

# elastic net with coarse age groups 
for (sex in c("female", "male")){
  for (ss in c("with", "without")){
    for (ag in coarse_age_groups) {
      #make sbatch
      rjobsh <- paste0("coarse_en_", sex, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
      
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                   "#SBATCH --time=12:00:00",
                   "#SBATCH --mem=60GB",
                   "#SBATCH --nodes=1",
                   "#SBATCH --cpus-per-task=1",
                   paste0("#SBATCH --job-name=", "coarse_en_", sex, "_", ag, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng"),
                   "#SBATCH --account=wang-krishnan",
                   "",
                   "cd $SLURM_SUBMIT_DIR",
                   "source /mnt/home/john3491/.bashrc",
                   "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                   "which python",
                   "",
                   paste0("python LR_elasticnet_3F_CV_for_one_age_group.py --sex=",
                          sex, " --age_group_type=coarse", " --pos_age_group=", ag,
                          " --data_type=rnaseq", " --std_scale=", ss, " --no-asinh"),
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

# l2 with fine age groups 
for (sex in c("female", "male")){
  for (ss in c("with", "without")){
    for (ag in fine_age_groups) {
      #make sbatch
      rjobsh <- paste0("l2_", sex, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
      
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                   "#SBATCH --time=4:00:00",
                   "#SBATCH --mem=60GB",
                   "#SBATCH --nodes=1",
                   "#SBATCH --cpus-per-task=1",
                   paste0("#SBATCH --job-name=", "l2_", sex, "_", ag, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng"),
                   "#SBATCH --account=wang-krishnan",
                   "",
                   "cd $SLURM_SUBMIT_DIR",
                   "source /mnt/home/john3491/.bashrc",
                   "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                   "which python",
                   "",
                   paste0("python LR_l2_3F_CV_for_one_age_group.py --sex=",
                          sex, " --age_group_type=fine", " --pos_age_group=", ag,
                          " --data_type=rnaseq", " --std_scale=", ss, " --no-asinh"),
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

# l2 with coarse age groups 
for (sex in c("female", "male")){
  for (ss in c("with", "without")){
    for (ag in coarse_age_groups) {
      #make sbatch
      rjobsh <- paste0("coarse_l2_", sex, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng", ".rjob.sh"); cat(rjobsh, "\n")
      
      rjobConn <- file(paste0(dirname,"/",rjobsh))
      writeLines(c("#!/bin/sh -login",
                   "#SBATCH --time=4:00:00",
                   "#SBATCH --mem=60GB",
                   "#SBATCH --nodes=1",
                   "#SBATCH --cpus-per-task=1",
                   paste0("#SBATCH --job-name=", "coarse_l2_", sex, "_", ag, "_", ag, "_rnaseq_no-asinh_", ss, "_stdsclng"),
                   "#SBATCH --account=wang-krishnan",
                   "",
                   "cd $SLURM_SUBMIT_DIR",
                   "source /mnt/home/john3491/.bashrc",
                   "export PATH=\"/mnt/home/john3491/anaconda3/bin:$PATH\"",
                   "which python",
                   "",
                   paste0("python LR_l2_3F_CV_for_one_age_group.py --sex=",
                          sex, " --age_group_type=coarse", " --pos_age_group=", ag,
                          " --data_type=rnaseq", " --std_scale=", ss, " --no-asinh"),
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
