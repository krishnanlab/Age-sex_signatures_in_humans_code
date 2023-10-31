# this script pulls metadata from all experimental directories and puts it together
# refine.bio RNA-seq compendia (HOMO_SAPIENS_2_1601762743.zip) unzips to one directory for each project, 
# with sample metadata and one file per sample in each directory
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)


project_dirs <- list.dirs("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio",
                          full.names = T, recursive = F)

#read in all metadata files (one per project), bind them together 
metadata_files <- paste0(project_dirs, "/metadata_", basename(project_dirs), ".tsv")
metadata <- lapply(metadata_files, read_delim, delim = "\t", col_names = T, col_types = cols(.default = "c"))
metadata <- bind_rows(metadata)

#write to file
metadata %>% 
  write_delim("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio/refine.bio_aggregated_metadata.tsv", 
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
