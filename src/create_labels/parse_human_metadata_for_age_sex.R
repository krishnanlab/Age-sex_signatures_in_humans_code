tic <- as.integer(as.POSIXct(Sys.time()))
library(tidyverse)

mdata <- read_delim("/mnt/research/compbio/krishnanlab/data/rnaseq/archs4/age_sex_project/filtered_samples_over_half_zero_exp/transformed_data/human_metadata.tsv",
                    delim = "\t", col_names = T)
#subset to useful columns
mdata <- mdata %>% 
  select(Sample_series_id, Sample_geo_accession, Sample_source_name_ch1,
         Sample_characteristics_ch1, Sample_description)
#make labels output
labels <- mdata %>% select(Sample_series_id, Sample_geo_accession)

###find males###
male_labels <- labels
#possible males
male_labels <- male_labels %>% 
  mutate(sn = ifelse(grepl("male", mdata$Sample_source_name_ch1, ignore.case = T), 
                     mdata$Sample_source_name_ch1, NA)) %>% 
  mutate(char = ifelse(grepl("male", mdata$Sample_characteristics_ch1, ignore.case = T), 
                     mdata$Sample_characteristics_ch1, NA)) %>% 
  mutate(desc = ifelse(grepl("male", mdata$Sample_description, ignore.case = T), 
                     mdata$Sample_description, NA))
#get rid of females
male_labels$sn <- ifelse(grepl("female", male_labels$sn, ignore.case = T), NA, male_labels$sn)
male_labels$char <- ifelse(grepl("female", male_labels$char, ignore.case = T), NA, male_labels$char)
male_labels$desc <- ifelse(grepl("female", male_labels$desc, ignore.case = T), NA, male_labels$desc)

#has male written in one of 3 cols
male_labels <- male_labels %>% 
  filter(!is.na(sn) | !is.na(char) | !is.na(desc))
#get one column T for male
male_labels <- male_labels %>% 
  select(Sample_series_id, Sample_geo_accession) %>% 
  mutate(male = T)


###find females###
female_labels <- labels
female_labels <- female_labels %>% 
  mutate(sn = ifelse(grepl("female", mdata$Sample_source_name_ch1, ignore.case = T), 
                     mdata$Sample_source_name_ch1, NA)) %>% 
  mutate(char = ifelse(grepl("female", mdata$Sample_characteristics_ch1, ignore.case = T), 
                       mdata$Sample_characteristics_ch1, NA)) %>% 
  mutate(desc = ifelse(grepl("female", mdata$Sample_description, ignore.case = T), 
                       mdata$Sample_description, NA))
#has female written in one of 3 cols
female_labels <- female_labels %>% 
  filter(!is.na(sn) | !is.na(char) | !is.na(desc))
#get one column T for female
female_labels <- female_labels %>% 
  select(Sample_series_id, Sample_geo_accession) %>% 
  mutate(female = T)

###find ages###
age_labels <- labels
#pull possible ages
age_labels$tmp <- mdata$Sample_characteristics_ch1
age_labels$tmp <- gsub("stage:", "", age_labels$tmp, ignore.case = T)
age_labels$tmp <- gsub("lineage:", "", age_labels$tmp, ignore.case = T)
age_labels$tmp <- gsub("passage:", "", age_labels$tmp, ignore.case = T)
#get "age:" without ".+age:"
age_labels <- age_labels %>% 
  mutate(char = ifelse(grepl("age:", tmp), tmp, NA)) %>% 
  select(-tmp) %>% 
  filter(!is.na(char))

#sep into fetus and after birth
fetus_age_labels <- age_labels %>% 
  filter(grepl("fetus", char, ignore.case = T) | 
           grepl("embryo", char, ignore.case = T) |
           grepl("gestation", char, ignore.case = T) |
           grepl("fetal", char, ignore.case = T) |
           grepl("Carnegie", char, ignore.case = T)) %>% 
  filter(!grepl("placenta", char, ignore.case = T)) %>% 
  filter(!grepl("mean", char, ignore.case = T)) %>% 
  select(-char) %>% 
  mutate(fetus = T)

post_birth_age_labels <- age_labels %>% 
  filter(!grepl("fetus", char, ignore.case = T) & 
           !grepl("embryo", char, ignore.case = T) &
           !grepl("gestation", char, ignore.case = T) &
           !grepl("fetal", char, ignore.case = T) &
           !grepl("Carnegie", char, ignore.case = T)) %>% 
  mutate(char = str_extract(char, "age: .{1,20}")) %>% 
  filter(!grepl("mean", char, ignore.case = T))
#remove "age: " from beginning of string
post_birth_age_labels$char <- sub("age: ", "", post_birth_age_labels$char)
###this filtering done with some manual pattern finding###
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!grepl("^NA", char)) %>% 
  filter(!grepl("^36-74", char)) %>% 
  filter(!grepl("^P", char, ignore.case = T)) %>% 
  filter(!grepl("^day", char, ignore.case = T)) %>% 
  filter(!grepl("^d", char, ignore.case = T)) %>% 
  filter(!grepl("^NES", char, ignore.case = T))

clean <- post_birth_age_labels %>% 
  filter(grepl("[0-9]+Xx.*$", char))
clean$char <- sub("Xx.*$", "", clean$char)
clean$char <- sub("4-6", "6", clean$char)
clean$char <- sub("3-5", "5", clean$char)
clean$char <- sub("16-19", "16", clean$char)
clean$char <- sub("15-30", "NA", clean$char)
clean$char <- sub("58-60", "58", clean$char)
clean$char <- sub("53-60", "53", clean$char)
clean$char <- sub("6-13", "NA", clean$char)
clean$char <- sub("6-8", "6", clean$char)
clean$char <- sub("8-10", "8", clean$char)
clean$char <- sub("8-12", "NA", clean$char)
clean$char <- sub("20-25", "25", clean$char)
clean$char <- sub("5-15", "NA", clean$char)
clean$char <- sub("23-28", "28", clean$char)
clean$char <- sub("1-3", "NA", clean$char)
clean$char <- sub("10-15", "15", clean$char)
clean$char <- sub("8-9", "8", clean$char)
clean$char <- sub("18-45", "NA", clean$char)
clean$char <- sub("40-45", "40", clean$char)
clean$char <- sub("30-33", "33", clean$char)
clean$char <- sub("5-10", "5", clean$char)
clean$char <- sub("27-34", "30", clean$char)
clean$char <- sub("Age ", "", clean$char)
clean$char <- sub("27/34", "30", clean$char)
clean$char <- sub("60 to 70", "60", clean$char)
clean$char <- sub("6 to 10", "8", clean$char)
clean$char <- sub("10 to 15", "12", clean$char)
clean$char <- sub("7 to 10", "8", clean$char)
clean$char <- sub("<20", "NA", clean$char)
clean$char <- sub("year ", "", clean$char)
clean$char <- sub("week[0-9]+", "NA", clean$char)

tmp <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% clean$Sample_geo_accession) %>% 
  filter(grepl("[0-9]+-[0-9]+", char)) %>% 
  mutate(char = str_extract(char, "[0-9]+-[0-9]+"))
tmp$char <- sub("22-36", "22", tmp$char)
tmp$char <- sub("64-88", "64", tmp$char)
tmp$char <- sub("22-40", "22", tmp$char)
tmp$char <- sub("21-30", "22", tmp$char)
tmp$char <- sub("30-40", "32", tmp$char)
tmp$char <- sub("27-45", "28", tmp$char)
tmp$char <- sub("21-34", "32", tmp$char)
tmp$char <- sub("12-18", "12", tmp$char)
#add to clean
clean <- bind_rows(clean, tmp)
#take out of unclean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% clean$Sample_geo_accession)
clean <- clean %>% 
  filter(!char == "NA")

tmp <- post_birth_age_labels %>% 
  filter(grepl("[0-9]+.*years", char, ignore.case = T)) %>% 
  mutate(char = tolower(char)) %>% 
  mutate(char = str_extract(char, "[0-9]+.*years"))
tmp$char <- sub(" years", "", tmp$char)
tmp$char <- sub("-years", "", tmp$char)
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean,tmp)

tmp <- post_birth_age_labels %>% 
  filter(!grepl("[^0-9]", char))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean,tmp)

tmp <- post_birth_age_labels %>% 
  filter(grepl("day", char, ignore.case = T)) %>% 
  mutate(char = ifelse(grepl("culture", char, ignore.case = T),
                       "NA", char)) %>% 
  mutate(char = ifelse(grepl("extract", char, ignore.case = T),
                       "NA", char)) %>% 
  mutate(char = tolower(char)) %>% 
  mutate(char = str_extract(char, "[0-9]+.*day"))
tmp$char <- sub(" day", "", tmp$char)
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
#clean to add to clean
tmp <- tmp %>% 
  filter(char != "NA") %>% 
  mutate(char = as.numeric(char)) %>% 
  mutate(char = round(char/365, digits = 2))
#make clean numeric
clean <- clean %>% 
  mutate(char = as.numeric(char))
clean <- bind_rows(clean, tmp)

tmp <- post_birth_age_labels %>% 
  filter(grepl("[0-9]+.*yrs", char, ignore.case = T)) %>% 
  mutate(char = str_extract(char, "[0-9]+")) %>%
  mutate(char = as.numeric(char))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean, tmp)

tmp <- post_birth_age_labels %>% 
  filter(grepl("[0-9]+.{0,2}y", char, ignore.case = T)) %>% 
  mutate(char = sub("<2", "1", char)) %>% 
  mutate(char = str_extract(char, "^[0-9]+")) %>%
  mutate(char = as.numeric(char))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean, tmp)

tmp <- post_birth_age_labels %>% 
  filter(grepl("w[e]*k", char, ignore.case = T)) %>% 
  mutate(char = str_extract(char, "^[0-9]+")) %>%
  mutate(char = as.numeric(char)) %>% 
  mutate(char = round(char/52, digits = 3))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean, tmp)
  
tmp <- post_birth_age_labels %>% 
  filter(grepl("^[0-9]+W", char)) %>% 
  mutate(char = str_extract(char, "^[0-9]+")) %>%
  mutate(char = as.numeric(char)) %>% 
  mutate(char = round(char/52, digits = 3))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean, tmp)

tmp <- post_birth_age_labels %>% 
  filter(grepl("month", char)) %>% 
  mutate(char = str_extract(char, "^[0-9]+")) %>%
  mutate(char = as.numeric(char)) %>% 
  mutate(char = round(char/12, digits = 3))
#take tmp out of post labels, add to clean
post_birth_age_labels <- post_birth_age_labels %>% 
  filter(!Sample_geo_accession %in% tmp$Sample_geo_accession)
clean <- bind_rows(clean, tmp)

colnames(clean) <- c("gse", "gsm", "age")
clean <- clean %>% 
  mutate(age_group = ifelse(age < 2, "infant", 
                            ifelse((age >= 2 & age < 10), "child",
                                   ifelse((age >= 10 & age < 20), "adolescent",
                                          ifelse((age >= 20 & age < 50), "adult",
                                                 ifelse((age >= 50 & age < 80), "olderadult", "elderly"))))))
fetus_age_labels <- fetus_age_labels %>% 
  select(-fetus) %>% 
  mutate(age = NA, age_group = "fetus")
colnames(fetus_age_labels) <- c("gse", "gsm", "age", "age_group")

#finish age labels
age_labels <- bind_rows(clean, fetus_age_labels)

#put everything together
colnames(male_labels) <- c("gse", "gsm", "sex")
colnames(female_labels) <- c("gse", "gsm", "sex")

male_labels <- male_labels %>% 
  mutate(sex = "male")
female_labels <- female_labels %>% 
  mutate(sex = "female")
#a few samples had both male and female in some combo of columns
bad_gsms <- male_labels %>% 
  filter(gsm %in% female_labels$gsm) %>% 
  pull(gsm)
  
male_labels <- male_labels %>% 
  filter(!gsm %in% bad_gsms)
female_labels <- female_labels %>% 
  filter(!gsm %in% bad_gsms)

sex_labels <- bind_rows(female_labels, male_labels)

all_labels <- full_join(sex_labels, age_labels, by = c("gse", "gsm"))
  
all_labels %>% 
  as.data.frame() %>% 
  write_delim("human_ARCHS4_text_parsed_gse-gsm-sex-age-age_group.tsv",
              delim = "\t", col_names = T)


#script time
toc <- as.integer(as.POSIXct(Sys.time()))
print(paste("The time it took to run this script in minutes was", (toc-tic)/60, sep = " "))



  
