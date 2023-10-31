# this script writes the ffq fetch commands for metadata
# for refine.bio samples
# ffq fetches SRA metadata: https://github.com/pachterlab/ffq
# cmd: ffq -t SRR -o srr_test.json SRR013580 SRR013560 SRR013569

tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

#read in refine.bio metadata to get run accessions
rb_run_accessions <- read_delim("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/refine.bio_aggregated_metadata.tsv",
                                delim = "\t", col_names = T, col_types = cols(.default = "c")) %>% 
  pull(refinebio_accession_code)

#sep DRR/ERR/SRR accessions (diff ffq commands)
drr_accessions <- rb_run_accessions[grepl("DRR", rb_run_accessions)]
err_accessions <- rb_run_accessions[grepl("ERR", rb_run_accessions)]
srr_accessions <- rb_run_accessions[grepl("SRR", rb_run_accessions)]

# split accessions into equal size groups (100 samples in each)
# output = list of length ceiling(total_runs/100)
drr_split <- split(drr_accessions, ceiling(seq_along(drr_accessions)/100))
err_split <- split(err_accessions, ceiling(seq_along(err_accessions)/100))
srr_split <- split(srr_accessions, ceiling(seq_along(srr_accessions)/100))

#collect ffq lines for drr runs in one vector for placement in sbatch later
drr_ffq_lines <- c("the first string")

for (run_list in drr_split){
  run_one <- head(run_list, n = 1)
  run_last <- tail(run_list, n = 1)
  ffq_lines <- paste0("ffq -t DRR -o drr_accessions_", run_one,
                      "_", run_last, ".json", " ", paste(run_list, collapse = " "))
  drr_ffq_lines <- c(drr_ffq_lines, ffq_lines)
}

drr_ffq_lines <- drr_ffq_lines[-1]

#collect ffq lines for err runs in one vector for placement in sbatch later
err_ffq_lines <- c("the first string")

for (run_list in err_split){
  run_one <- head(run_list, n = 1)
  run_last <- tail(run_list, n = 1)
  ffq_lines <- paste0("ffq -t ERR -o err_accessions_", run_one,
                      "_", run_last, ".json", " ", paste(run_list, collapse = " "))
  err_ffq_lines <- c(err_ffq_lines, ffq_lines)
}

err_ffq_lines <- err_ffq_lines[-1]

#collect ffq lines for srr runs in one vector for placement in sbatch later
srr_ffq_lines <- c("the first string")

for (run_list in srr_split){
  run_one <- head(run_list, n = 1)
  run_last <- tail(run_list, n = 1)
  ffq_lines <- paste0("ffq -t SRR -o srr_accessions_", run_one,
                      "_", run_last, ".json", " ", paste(run_list, collapse = " "))
  srr_ffq_lines <- c(srr_ffq_lines, ffq_lines)
}

srr_ffq_lines <- srr_ffq_lines[-1]

#create directory for sbatches
dirname <- "./ffq_sbatches"
if(!dir.exists(dirname)) {
  dir.create(dirname)
}
#write drr lines to sbatch
drr_file <- file("./ffq_sbatches/drr_ffq_commands.sbatch")
writeLines(c("#!/bin/sh -login",
             "#SBATCH --time=14:00:00",
             "#SBATCH --mem=10GB",
             "#SBATCH --nodes=1",
             "#SBATCH --cpus-per-task=1",
             "#SBATCH --job-name=drr_ffq",
             "#SBATCH --account=wang-krishnan",
             "",
             "cd $SLURM_SUBMIT_DIR",
             "",
             drr_ffq_lines), 
           con = drr_file)
close(drr_file)


#ten hours didn't always cut it for 25 ffq commands 
#of 100 samples each, so changed script to write 14 hour jobs
#split err lines into groups of 25
err_line_split <- split(err_ffq_lines, ceiling(seq_along(err_ffq_lines)/25))
#write err lines to several sbatches 
for (index in seq_along(err_line_split)){
  err_file <- file(paste0("./ffq_sbatches/err_ffq_commands_file", index, ".sbatch"))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=14:00:00",
               "#SBATCH --mem=10GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=err_", index),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "",
               err_line_split[[index]]), 
             con = err_file)
  close(err_file)
}

#split srr lines into groups of 25
srr_line_split <- split(srr_ffq_lines, ceiling(seq_along(srr_ffq_lines)/25))
#write srr lines to several sbatches 
for (index in seq_along(srr_line_split)){
  srr_file <- file(paste0("./ffq_sbatches/srr_ffq_commands_file", index, ".sbatch"))
  writeLines(c("#!/bin/sh -login",
               "#SBATCH --time=14:00:00",
               "#SBATCH --mem=10GB",
               "#SBATCH --nodes=1",
               "#SBATCH --cpus-per-task=1",
               paste0("#SBATCH --job-name=srr_", index),
               "#SBATCH --account=wang-krishnan",
               "",
               "cd $SLURM_SUBMIT_DIR",
               "",
               srr_line_split[[index]]), 
              con = srr_file)
  close(srr_file)
}

#submit all sbatches
sbatches <- list.files("./ffq_sbatches", full.names = F)
for (sbatch in sbatches){
  system(paste0("sbatch ", paste0("./ffq_sbatches/", sbatch)))
  
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