tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

src_dir <- "/mnt/home/john3491/projects/age-sex-prediction/src/age_prediction_by_sex/"

# make a directory to put sbatches
dirname <- paste0("./sbatches_term_overlap")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

data_path <- "~/projects/age-sex-prediction/data/gene_enrichment_for_sex-biased_analyses/"
file_names <- c("RNAseq_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv",
                "microarray_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv",
                "sagd_gene-group_file.tsv", 
                "gtex-genes_gene-group_file.tsv",
                "guo_gene-group_file.tsv")

for (file_name in file_names) {
        #make sbatch
        rjobsh <- paste0("rnaseq_overlap", file_name, ".rjob.sh"); cat(rjobsh, "\n")
        
        rjobConn <- file(paste0(dirname,"/",rjobsh))
        writeLines(c("#!/bin/sh -login",
                     "#SBATCH --time=2:00:00",
                     "#SBATCH --mem=30GB",
                     "#SBATCH --nodes=1",
                     "#SBATCH --cpus-per-task=1",
                     paste0("#SBATCH --job-name=rnaseq_overlap_", file_name),
                     "#SBATCH --account=wang-krishnan",
                     "",
                     "cd $SLURM_SUBMIT_DIR",
                     "module purge",
                     "module load GCC/8.3.0 OpenMPI/3.1.4",
                     "module load R/4.0.2",
                     "",
                     paste("Rscript term_enrichment_between_multiple_sets.R",
                           paste0(data_path, "RNAseq_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv"), 
                           paste0(data_path, file_name), sep = " "),
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

for (file_name in file_names) {
  #make sbatch
  rjobsh <- paste0("microarray_overlap", file_name, ".rjob.sh"); cat(rjobsh, "\n")
  
  rjobConn <- file(paste0(dirname,"/",rjobsh))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=2:00:00",
               "#SBATCH --mem=30GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=microarray_overlap_", file_name),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "module purge",
               "module load GCC/8.3.0 OpenMPI/3.1.4",
               "module load R/4.0.2",
               "",
               paste("Rscript term_enrichment_between_multiple_sets.R",
                     paste0(data_path, "microarray_age_group_elasticnet_models_genes_over_2sd_from_mean.tsv"), 
                     paste0(data_path, file_name), sep = " "),
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

