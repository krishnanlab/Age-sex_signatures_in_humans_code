tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)


# make a directory to put sbatches
dirname <- paste0("./sbatches_sex-strat_gtex_med-tpm")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

# tissue with female samples
fem_tis <- c("Blood", "Brain", "Adipose", "Muscle", "BV",
             "Heart", "Ovary", "Uterus", "Vagina", "Breast", "Skin", "Salivary",
             "Lung", "Thyroid", "Adrenal", "Spleen", "Pancreas", "Small",
             "Colon", "Esophagus", "Stomach", "Pituitary", "Nerve", "Liver",
             "Kidney", "Cervix", "Fallopian Tube", "Bladder", "Bone")
# tissue with male samples
mal_tis <- c("Blood", "Adrenal", "Thyroid", "Lung", "Spleen", "Pancreas",
             "Esophagus", "Stomach", "Adipose", "Skin", "Colon", "Small",
             "Prostate", "Testis", "Muscle", "Nerve", "Brain", "BV",
             "Heart", "Salivary", "Pituitary", "Breast", "Liver", "Kidney",
             "Bladder")

arg_list <- c(paste("female", fem_tis, sep = " "),
              paste("male", mal_tis, sep = " "))


for (arg_set in arg_list) {
  #make sbatch
  rjobsh <- paste0("sex-strat_med-tpm_", gsub(" ", "_", arg_set), ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=4:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=sex-strat_med-tpm_", gsub(" ", "_", arg_set)),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "module purge",
               "module load GCC/8.3.0 OpenMPI/3.1.4",
               "module load R/4.0.2",
               "",
               paste("Rscript multi-job_sex-strat_gtex_median_gene_value_per_tissue.R",
                     arg_set, sep = " "),
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
