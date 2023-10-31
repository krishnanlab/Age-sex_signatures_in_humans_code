tic <- as.integer(as.POSIXct(Sys.time()))

library(tidyverse)
library(jsonlite)

# list files in ffq dir
file_list <- list.files("./ffq", full.names = T, pattern = "*json")
# read each into a list (which makes files a list of lists)
files <- lapply(file_list, read_json)

# read in rs labels
rs_labels <- read_delim("~/projects/age-sex-prediction/data/labels/metasra-refine.bio_sample-run-experiment-sex-age-tissue-cell_line-disease.tsv",
                        delim = "\t", col_names = T)

# create empty list to add to
sample_attributes <- list()
# loop over list of lists to pull sample attributes out
for (file in files){
  tmp <- lapply(file, pluck, "sample", "attributes")
  sample_attributes <- c(sample_attributes, tmp)
}

# stack list/bind rows into dataframe
sample_attributes <- lapply(sample_attributes, stack)  

sample_attributes <- bind_rows(sample_attributes, .id = "run")

sample_attributes <- as_tibble(sample_attributes)

# remove ENA stuff, won't help in checking age/sex labels
# lots of other random stuff would take too long to remove
cell_line_runs <- sample_attributes %>% 
  filter(grepl("cell", ind)) %>% 
  filter(grepl("line", ind)) %>% 
  filter(values != "not applicable") %>% 
  filter(values != "Left frontal") %>% 
  filter(values != "Left occipital") %>%
  filter(values != "Right frontal") %>%
  filter(values != "Left parietotemporal") %>%
  filter(values != "Left parietal") %>%
  filter(values != "Left temporal") %>%
  filter(values != "Right parietal") %>%
  filter(values != "Right temporo-occipital") %>%
  filter(values != "Right frontotemporal") %>%
  filter(values != "Multifocal (left frontal/temporal/parietal; right occipital)") %>%
  filter(values != "Right temporal") %>%
  filter(values != "Right frontal/basal ganglia") %>%
  filter(values != "Right parietal; right temporal") %>% 
  filter(values != "primary skeletal muscle") %>% 
  filter(values != "primary ovary") %>%
  filter(values != "primary breast") %>%
  filter(values != "primary naive B cells") %>%
  filter(values != "primary testis") %>% 
  filter(values != "primary tissue") %>% 
  pull(run)

cl_runs <- sample_attributes %>% 
  filter(grepl("cell", values, ignore.case = T)) %>% 
  filter(grepl("line", values, ignore.case = T)) %>% 
  pull(run)

cc_runs <- sample_attributes %>% 
  filter(grepl("culture", ind, ignore.case = T) | grepl("culture", ind, ignore.case = T)) %>% 
  filter(!grepl("not", values, ignore.case=T)) %>% 
  filter(!grepl("none", values, ignore.case=T)) %>% 
  filter(values != "No") %>% 
  pull(run)

pass_runs <- sample_attributes %>% 
  filter(grepl("passage", ind)) %>% 
  filter(!grepl("not applicable", values, ignore.case = T)) %>% 
  pull(run)

sc_runs <- sample_attributes %>% 
  filter(ind == "single cell well quality") %>% 
  pull(run)

mouse_runs <- sample_attributes %>% 
  filter(grepl("mouse", ind)) %>% 
  filter(values != "NA") %>% 
  pull(run)

tf_runs <- sample_attributes %>% 
  filter(grepl("transfect", ind)) %>% 
  filter(values != "untransfected") %>% 
  filter(values != "Untransfected") %>% 
  filter(values != "none") %>% 
  pull(run)

mixed_runs <- sample_attributes %>% 
  filter(grepl("mixed sample", values, ignore.case=T)) %>% 
  pull(run)

graft_runs <- sample_attributes %>% 
  filter(grepl("graft", values)) %>% 
  pull(run)

g_runs <- sample_attributes %>% 
  filter(grepl("graft", ind)) %>% 
  filter(values != "false") %>% 
  filter(values != "not applicable") %>% 
  filter(values != "No")

tp_runs <- sample_attributes %>% 
  filter(grepl("transplant", values)) %>% 
  pull(run)

t_runs <- sample_attributes %>% 
  filter(grepl("transplant", ind)) %>% 
  pull(run)

unacceptable_runs <- unique(c(cell_line_runs, cl_runs, cc_runs, pass_runs, sc_runs, mouse_runs, 
                              tf_runs, graft_runs, g_runs, tp_runs, mixed_runs))

sample_attributes <- sample_attributes %>% 
  filter(!grepl("ENA-*", ind)) %>% 
  filter(ind != "sample_name") %>% 
  filter(ind != "bioproject_id") %>% 
  filter(ind != "isolate") %>% 
  filter(ind != "ID") %>% 
  filter(ind != "organism") %>%
  filter(ind != "BioSampleModel") %>% 
  filter(ind != "biomaterial_provider") %>% 
  filter(!(run %in% unacceptable_runs))

sample_attributes <- sample_attributes %>% 
  group_by(run) %>% 
  summarise(values = paste(values, collapse = "|"), ind = paste(ind, collapse = "|"))

# start looking for fetus labels
# metaSRA did not have any samples labeled for fetus
gest <- sample_attributes %>% 
  filter(grepl("gestation", values))
fetus <- sample_attributes %>% 
  filter(grepl("fetus", values))
fetal <- sample_attributes %>% 
  filter(grepl("fetal", values))

poss_fetus <- bind_rows(gest, fetus)
poss_fetus <- bind_rows(poss_fetus, fetal)

gest <- sample_attributes %>% 
  filter(grepl("gestation", ind))
fetus <- sample_attributes %>% 
  filter(grepl("fetus", ind))
fetal <- sample_attributes %>% 
  filter(grepl("fetal", ind))

poss_fetus <- bind_rows(poss_fetus, gest)
poss_fetus <- bind_rows(poss_fetus, fetus)
poss_fetus <- bind_rows(poss_fetus, fetal)

poss_fetus <- poss_fetus[!duplicated(poss_fetus),]

# get rid of the 13 runs that were pregnant women
poss_fetus <- poss_fetus %>% 
  filter(!grepl("pregnant", values)) %>% 
  filter(!grepl("gestation 37+2 week", values, fixed = T))

fetus_labels <- poss_fetus %>% 
  separate_rows(values, ind, sep = "\\|")

age_fields <- c("weeks of gestation", "gestational week",
                "age of the fetus", "development stage", 
                "Stage", "gestational age", "dev_stage", "age",
                "developmental stage")
age_labels <- fetus_labels %>% 
  filter(ind %in% age_fields) %>% 
  filter(!grepl("neonate", values)) %>% 
  filter(!grepl("cell", values)) %>% 
  filter(!(ind == "gestational age" & values == "other")) %>% 
  mutate(age_group = "fetus") %>%  
  mutate(age = values)

age_labels$age <- gsub("54-57 days post conception (gestation)", "8", age_labels$age, fixed = T)
age_labels$age <- gsub("fetal (12-18 gestational weeks)", "12-18", age_labels$age, fixed = T)
age_labels$age <- gsub("16 weeks 6 days", "17", age_labels$age)
age_labels$age <- gsub("15 weeks 6 days", "16", age_labels$age)
age_labels$age <- gsub("91 days", "13", age_labels$age)
age_labels$age <- gsub("103 days", "15", age_labels$age)
age_labels$age <- gsub("110 days", "16", age_labels$age)
age_labels$age <- gsub("114 days", "16", age_labels$age)
age_labels$age <- gsub("6w 5d", "7", age_labels$age)
age_labels$age <- gsub("5 and 6 week gestation", "5.5", age_labels$age)
age_labels$age <- gsub("prenatal", "", age_labels$age)
age_labels$age <- gsub("25 3/7 weeks", "25", age_labels$age, fixed = T)
age_labels$age <- gsub("26 2/7 weeks", "26", age_labels$age, fixed = T)
age_labels$age <- gsub("27 2/7 weeks", "27", age_labels$age, fixed = T)
age_labels$age <- gsub("24 6/7 weeks", "25", age_labels$age, fixed = T)
age_labels$age <- gsub("25 4/7 weeks", "26", age_labels$age, fixed = T)
age_labels$age <- gsub("29 6/7 weeks", "30", age_labels$age, fixed = T)
age_labels$age <- gsub("29 2/7 weeks", "29", age_labels$age, fixed = T)
age_labels$age <- gsub("33 5/7 weeks", "34", age_labels$age, fixed = T)
age_labels$age <- gsub("26 6/7 weeks", "27", age_labels$age, fixed = T)
age_labels$age <- gsub("26 5/7 weeks", "27", age_labels$age, fixed = T)
age_labels$age <- gsub("33 6/7 weeks", "34", age_labels$age, fixed = T)
age_labels$age <- gsub("16 weeks 3 days", "17", age_labels$age)
age_labels$age <- gsub("16 weeks 4 days", "17", age_labels$age)
age_labels$age <- gsub("17 weeks 2 days", "18", age_labels$age)
age_labels$age <- gsub("post-conception", "", age_labels$age, fixed = T)
age_labels$age <- gsub("gestational", "", age_labels$age, ignore.case = T)
age_labels$age <- gsub("gestation", "", age_labels$age, ignore.case = T)
age_labels$age <- gsub("weeks", "", age_labels$age, ignore.case = T)
age_labels$age <- gsub("week", "", age_labels$age, ignore.case = T)
age_labels$age <- gsub("(pcw)", "", age_labels$age, fixed = T)
age_labels$age <- gsub(" to ", "-", age_labels$age, fixed = T)
age_labels$age <- gsub("postconceptional", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("fetal", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("fetus", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("age", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("of", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("after", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("liver", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("w", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub("k", "", age_labels$age,ignore.case = T)
age_labels$age <- gsub(",", "", age_labels$age)
age_labels$age <- gsub("th", "", age_labels$age)
age_labels$age <- gsub("s", "", age_labels$age)
age_labels$age <- trimws(age_labels$age, which = "both")
age_labels$age <- paste0(age_labels$age, "wk")

age_labels <- age_labels %>% 
  select(-ind, -values) %>% 
  filter(!(run %in% unacceptable_runs)) %>% 
  filter(age != "wk")

age_labels <- age_labels[!duplicated(age_labels),]

sex_labels <- fetus_labels %>% 
  filter(grepl("sex", ind, ignore.case = T) | grepl("gender", ind, ignore.case = T)) %>% 
  rename(sex = values) %>% 
  select(-ind) %>% 
  filter(!(run %in% unacceptable_runs))

sex_labels$sex <- gsub("N/A", NA, sex_labels$sex)
sex_labels$sex <- gsub("unknown", NA, sex_labels$sex, ignore.case = T)
sex_labels$sex <- gsub("not determined", NA, sex_labels$sex)
sex_labels$sex <- gsub("pooled male and female", "mixed", sex_labels$sex)
sex_labels$sex <- gsub("Male", "male", sex_labels$sex)
sex_labels$sex <- gsub("Female", "female", sex_labels$sex)

tissue_labels <- fetus_labels %>% 
  filter(grepl("tissue", ind, ignore.case = T) | grepl("source", ind, ignore.case = T) | grepl("region", ind, ignore.case = T)) %>% 
  filter(ind != "tissue archive method") %>% 
  rename(tissue = values) %>% 
  select(-ind) %>% 
  filter(!(run %in% unacceptable_runs)) %>% 
  filter(tissue != "normal fetus") %>% 
  filter(!(grepl("placenta", tissue, ignore.case = T))) %>% 
  filter(!(grepl("weeks fetal", tissue))) %>% 
  filter(!(grepl("post-gestation", tissue))) %>% 
  filter(tissue != "preculture") %>% 
  filter(tissue != "fetal human cortex 1") %>% 
  filter(tissue != "fetal human cortex 2") %>% 
  filter(!grepl("explants", tissue)) %>% 
  filter(!grepl("euploid", tissue)) %>% 
  filter(!grepl("twin", tissue)) %>% 
  filter(tissue != "human fetus") %>% 
  filter(tissue != "human embryonic retinal cells") %>% 
  filter(!grepl("_embryo", tissue)) %>% 
  filter(tissue != "Human fetal pancreata") %>% 
  filter(tissue != "embryonic renal corpuscles") %>% 
  filter(!grepl("gestational", tissue))
  

tissue_labels$tissue <- tolower(tissue_labels$tissue)
tissue_labels$tissue <- gsub("fetal ", "", tissue_labels$tissue)
tissue_labels$tissue <- gsub(", fetal", "", tissue_labels$tissue)

tissue_labels <- tissue_labels[!duplicated(tissue_labels),]

tissue_labels <- tissue_labels %>% 
  group_by(run) %>% 
  summarise(tissue = paste(tissue, collapse = ", "))

cell_labels <- fetus_labels %>% 
  filter(ind %in% c("cell type", "cell_type", "cell subtype", "cell_subtype")) %>% 
  rename(cells = values) %>% 
  group_by(run) %>% 
  summarise(cells = paste(cells, collapse = ", "))

disease_labels <- fetus_labels %>% 
  filter(grepl("disease", ind, ignore.case = T) | grepl("health", ind, ignore.case = T) | grepl("pathology", ind, ignore.case = T)) %>% 
  rename(disease = values) %>% 
  select(-ind) %>% 
  filter(!(run %in% unacceptable_runs))

age_labels <- left_join(age_labels, sex_labels, by = "run")
age_labels <- left_join(age_labels, tissue_labels, by = "run")
age_labels <- left_join(age_labels, cell_labels, by = "run")
age_labels <- left_join(age_labels, disease_labels, by = "run")

# get sample and experiment IDs for the possible fetus-labeled runs
fetus_runs <- unique(age_labels$run)
# get names of runs contained in each file
file_runs <- lapply(files, names)
names(file_runs) <- basename(file_list)

check_runs <- function(list){
  if (sum(fetus_runs %in% list) > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# return T/F for each file 
# (does it contain a possible fetus-labeled run?)
fetus_files <- lapply(file_runs, check_runs)
# turn into df
fetus_files <- stack(fetus_files)
fetus_files <- fetus_files %>% 
  filter(values == TRUE) %>% 
  pull(ind) %>% 
  as.character()
fetus_files <- paste0("./ffq/", fetus_files)
# read each into a list (which makes files a list of lists)
ffiles <- lapply(fetus_files, read_json)

# sample IDs
fetus_run_samples <- list()
for (file in ffiles){
  tmp <- lapply(file, pluck, "sample", "accession")
  fetus_run_samples <- c(fetus_run_samples, tmp)
}
fetus_run_samples <- stack(fetus_run_samples)
colnames(fetus_run_samples) <- c("sample", "run")
# add sample IDs to run labels
age_labels <- left_join(age_labels, fetus_run_samples, by = "run")

# experiment IDs
fetus_run_experiments <- list()
for (file in ffiles){
  tmp <- lapply(file, pluck, "study", "accession")
  fetus_run_experiments <- c(fetus_run_experiments, tmp)
}
fetus_run_experiments <- stack(fetus_run_experiments)
colnames(fetus_run_experiments) <- c("experiment", "run")
# add experiment IDs to run labels
age_labels <- left_join(age_labels, fetus_run_experiments, by = "run")

# change colnames of rs_labels
colnames(rs_labels) <- c("sample", "run", "experiment", 
                         "sex", "age", "age_group", 
                         "tissue", "cells", "disease")
# change order of age_labels (fetus labels) to match rs_labels
age_labels <- age_labels %>% select(sample, run, experiment, sex, age, age_group, tissue, cells, disease)

# join together
# make age character for rs_labels, because fetus labels have "wk" attached
rs_labels <- rs_labels %>% 
  arrange(age)
rs_labels <- rs_labels %>% 
  mutate(age = as.character(age)) %>% 
  filter(!(run %in% unacceptable_runs))

rs_labels <- bind_rows(age_labels, rs_labels)

# remove the few fetal samples that were tagged by
# metaSRA and labeled as age = 0 --> infant 
rs_labels <- rs_labels %>%
  filter(!(run == "SRR1044510" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044511" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044512" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044513" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044514" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044516" & age_group == "infant")) %>% 
  filter(!(run == "SRR1044517" & age_group == "infant")) 

# get rid of bad samples from metasra labeled samples
rs_labels %>% 
  write_delim("~/projects/age-sex-prediction/data/labels/refine.bio_sample-run-experiment-sex-age-tissue-cell_line-disease.tsv",
              delim = "\t", col_names = T)

# subset sample attributes to metaSRA labeled samples
sample_attributes <- sample_attributes %>% 
  filter(run %in% rs_labels$run)

# order runs in same order as rs_labels
reorder_idx <- match(rs_labels$run, sample_attributes$run) 

sample_attributes <- sample_attributes[reorder_idx,]

sample_attributes %>% 
  write_delim("~/projects/age-sex-prediction/data/labels/ffq_sample_attributes_to_check_metasra-refine.bio_labels.tsv",
              delim = "\t", col_names = T)

#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))
