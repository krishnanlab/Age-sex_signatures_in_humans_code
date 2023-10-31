# this script pulls metadata from all experimental directories and puts it together
# refine.bio RNA-seq compendia (HOMO_SAPIENS_2_1601762743.zip) unzips to one directory for each project, 
# with sample metadata and one file per sample in each directory
# this script also matches run IDs to sample IDs using
# the file written from the ffq fetches
# so metaSRA labels can be matched to refine.bio samples (all SRA)
# matched metaSRa labels will be written to this metadata file
tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)
library(jsonlite)

project_dirs <- list.dirs("/mnt/research/compbio/krishnanlab/data/rnaseq/refine.bio",
                      full.names = T, recursive = F)

# read in all metadata files (one per project), bind them together 
metadata_files <- paste0(project_dirs, "/metadata_", basename(project_dirs), ".tsv")
metadata <- lapply(metadata_files, read_delim, delim = "\t", col_names = T, col_types = cols(.default = "c"))
metadata <- bind_rows(metadata)
# rename column
metadata <- metadata %>% rename(run = refinebio_accession_code)

# read in sample/run map from ffq
id_conversion <- read_delim("~/projects/age-sex-prediction/data/refine.bio/refine.bio_runs_matched_to_samples.tsv",
                            delim = "\t", col_names = T)

# read in metasra labels json file
metasra <- fromJSON("/mnt/research/compbio/krishnanlab/data/metaSRA/metasra.v1-8.json")
# subset metasra file to sample ids actually in refine.bio
metasra <- metasra[names(metasra) %in% unique(id_conversion$sample)]

######get age labels############
# use purrr to pluck the real-value properties dataframes from each sample in json (result is named list)
real_values <- lapply(metasra, pluck, "real-value properties")
# bind together into one dataframe using names to id sample
real_values <- bind_rows(real_values, .id = "sample") %>% as_tibble()
# get labels for age (EFO:0000246)
real_values <- real_values %>% 
  filter(property_id == "EFO:0000246") %>% 
  select(-property_id)
#rename columns
colnames(real_values) <- c("sample", "metasra_age", "metasra_age_unit")
# replace ontology terms with words
# real_values$metasra_age_unit <- gsub("missing", NA, real_values$metasra_age_unit)
real_values$metasra_age_unit <- gsub("UO:0000036", "year", real_values$metasra_age_unit)
real_values$metasra_age_unit <- gsub("UO:0000035", "month", real_values$metasra_age_unit)
real_values$metasra_age_unit <- gsub("UO:0000034", "week", real_values$metasra_age_unit)
real_values$metasra_age_unit <- gsub("UO:0000033", "day", real_values$metasra_age_unit)
real_values$metasra_age_unit <- gsub("UO:0000032", "hour", real_values$metasra_age_unit)
# make ages in years
real_values <- real_values %>% 
  filter(!grepl("hour", metasra_age_unit)) %>% 
  mutate(metasra_age = ifelse(metasra_age_unit == "day", (metasra_age/365),
                              ifelse(metasra_age_unit == "week", (metasra_age/52),
                              ifelse(metasra_age_unit == "month", (metasra_age/12),
                              metasra_age)))) %>% 
  select(-metasra_age_unit)
# get rid of samples that were tagged with age more than once
multi_tagged <- real_values$sample[duplicated(real_values$sample)]
multi_tagged <- unique(multi_tagged)
real_values <- real_values %>% 
  filter(!sample %in% multi_tagged)
# make age_group column
real_values <- real_values %>% 
  mutate(metasra_age_group = ifelse(metasra_age < 0, "fetus",
         ifelse((metasra_age >= 0 & metasra_age <= 2), "infant",
                ifelse((metasra_age > 2 & metasra_age <= 10), "child",
                       ifelse((metasra_age > 10 & metasra_age <= 20), "adolescent",
                              ifelse((metasra_age > 20 & metasra_age <= 50), "adult",
                                     ifelse((metasra_age > 50 & metasra_age <= 80), "olderadult",
                                            ifelse(is.na(metasra_age), NA, "elderly"))))))))

########get sex and disease labels############
# this gets ont_values into tidy dataframe with sample | term columns
ont_values <- lapply(metasra, pluck, "mapped ontology terms")
ont_values <- plyr::ldply(ont_values, rbind) %>% as_tibble()
# rename, put terms into one column, separate into rows to make tidy, get rid of induced "NA"
# represented ontologies: DOID, UBERON, EFO (experimental factor), CL, and CVCL (cellosaurus)
ont_values <- ont_values %>% 
  rename(sample = .id) %>% 
  unite("metasra_terms", 2:31, sep = ", ") %>% 
  separate_rows(metasra_terms, sep = ", ") %>% 
  filter(metasra_terms != "NA")

# pull out sex labels
sex_labels <- ont_values %>% 
  filter(metasra_terms %in% c("UBERON:0003100","UBERON:0003101")) %>% 
  rename(metasra_sex = metasra_terms)
# match format of microarray labels
sex_labels$metasra_sex <- gsub("UBERON:0003100", "female", sex_labels$metasra_sex)
sex_labels$metasra_sex <- gsub("UBERON:0003101", "male", sex_labels$metasra_sex)
# get rid of samples that were tagged with sex twice
twice_tagged <- sex_labels$sample[duplicated(sex_labels$sample)]
twice_tagged <- unique(twice_tagged)
sex_labels <- sex_labels %>% 
  filter(!sample %in% twice_tagged)

# pull out disease labels
disease_labels <- ont_values %>% 
  filter(grepl("DOID", metasra_terms)) %>% 
  group_by(sample) %>% 
  summarise(metasra_disease_terms = paste(metasra_terms, collapse = ", "))

# pull out tissue labels 
# getting rid of sex since I know the terms,
# but remaining terms aren't necessarily all tissue terms
tissue_labels <- ont_values %>% 
  filter(grepl("UBERON", metasra_terms)) %>% 
  filter(!metasra_terms %in% c("UBERON:0003100","UBERON:0003101")) %>% 
  group_by(sample) %>% 
  summarise(metasra_tissue_terms = paste(metasra_terms, collapse = ", "))

# cell line labels
cell_labels <- ont_values %>% 
  filter(grepl("CL:", metasra_terms)) %>% 
  group_by(sample) %>% 
  summarise(metasra_cell_line_terms = paste(metasra_terms, collapse = ", "))

# add age/sex/ontology labels to id_conversion
id_conversion <- left_join(id_conversion, sex_labels, by = "sample")
id_conversion <- left_join(id_conversion, real_values, by = "sample")
id_conversion <- left_join(id_conversion, tissue_labels, by = "sample")
id_conversion <- left_join(id_conversion, cell_labels, by = "sample")
id_conversion <- left_join(id_conversion, disease_labels, by = "sample")

# add to metadata
metadata <- left_join(metadata, id_conversion, by = "run")

# write full data to file
metadata %>% 
  write_delim("~/projects/age-sex-prediction/data/refine.bio/full_refine.bio-metasra_aggregated_metadata.tsv", 
              delim = "\t", col_names = T)

# write age+sex labeled samples only 
# with only necessary columns
metadata %>% 
  filter(!is.na(metasra_age)) %>% 
  filter(!is.na(metasra_sex)) %>% 
  select(sample, run, experiment_accession, starts_with("metasra")) %>% 
  write_delim("~/projects/age-sex-prediction/data/labels/metasra-refine.bio_sample-run-experiment-sex-age-tissue-cell_line-disease.tsv",
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))

