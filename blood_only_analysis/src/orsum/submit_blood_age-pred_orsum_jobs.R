tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

term_sets <- c('gobp', 'gwas', 'mp', 'mondo')
gmt_dir <- "../../../zscore_gene_enrichment_analysis/data/orsum_gene_set_gmts"
gmt_files <- list.files(gmt_dir, full.names = T)
significance_dir <- "../../data/orsum_blood_age_prediction_significance_files/"
significance_files <- list.files(significance_dir, full.names = T)
output_dir <- "../../results/orsum/age-pred/"

# make a directory to put sbatches
dirname <- "./sbatches_age-pred_orsum_jobs"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

for (term_set in term_sets) {
  #make sbatch
  rjobsh <- paste0("orsum", term_set, ".rjob.sh"); cat(rjobsh, "\n")
  
  #get correct gmt
  gmt <- gmt_files[grepl(term_set, gmt_files, ignore.case = T)]
  
  sig_files <- significance_files[grepl(term_set, significance_files)]
  sig_files_com <- sQuote(sig_files, q = FALSE)
  sig_files_com <- paste(sig_files_com, collapse = " ")
  # write file
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=4:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=orsum_age_", term_set),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "source /mnt/home/john3491/.bashrc",
               "export PATH=/mnt/home/john3491/anaconda3/bin:$PATH",
               "which python",
               "",
               "conda activate orsum",
               paste0("orsum.py --gmt \'", gmt, "\' --files ",
                      sig_files_com,
                      " --outputFolder \'", output_dir,
                      term_set, "\'"),
               "",
               "echo \"finished job\""), rjobConn)
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

