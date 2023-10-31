#process json files from ffq fetch commands
#want sample | run conversion table so the
#refine.bio runs can be matched to 
#metaSRA sample labels
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(jsonlite)
#list files in ffq dir
file_list <- list.files("./ffq", full.names = T, pattern = "*json")
#read each into a list (which makes files a list of lists)
files <- lapply(file_list, read_json)

#create empty list to add to
run_samples <- list()
#loop over list of lists to pull sample accession out
for (file in files){
  tmp <- lapply(file, pluck, "sample", "accession")
  run_samples <- c(run_samples, tmp)
}

#stack named list into dataframe
run_samples <- stack(run_samples)
colnames(run_samples) <- c("sample", "run")

#write to file
run_samples %>% 
  write_delim("~/projects/age-sex-prediction/data/refine.bio/refine.bio_runs_matched_to_samples.tsv",
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
