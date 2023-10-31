#combine samples_not_sexed_by_metadata_parsing.tsv [output from predict_sex_for_ARCHS4_samples_with_extracted_age_info.py]
#and human_ARCHS4_text_parsed_gse-gsm-sex-age-age_group.tsv [output from parse_human_metadata_for_age_sex.R]

library(tidyverse)
path <- "/mnt/research/compbio/krishnanlab/data/rnaseq/archs4/age_sex_project/filtered_samples_over_half_zero_exp/transformed_data/"


unsexed <- read_delim(paste0(path, "samples_not_sexed_by_metadata_parsing.tsv"), delim = "\t", col_names = T)
#rename because it was not capitalized in original file
unsexed <- rename(unsexed, gsm = GSM)
#rename to more accurate - under cut on any gene == 'female' by that gene's cutoff,
#over cut on any gene == 'male' by that gene's cutoff
#therefore 0 == male by all gene cuts, 7 == 'female' by all gene cuts
unsexed <- rename(unsexed, number_genes_under_cut = sex)
all <- read_delim(paste0(path, "human_ARCHS4_text_parsed_gse-gsm-sex-age-age_group.tsv"), delim = "\t", col_names = T)

#only going to declare 6/7s as female, 0/1s as male
unsexed <- unsexed %>% 
  mutate(sex = ifelse(number_genes_under_cut > 5, 'female', 
                      ifelse(number_genes_under_cut < 2, 'male', NA)))

unsexed <- unsexed %>% select(gse, gsm, sex, age, age_group)

#just drop unsexed from all and then row bind
all <- all %>% filter(!gsm %in% unsexed$gsm)
#get rid of still NAs for sex
unsexed <- unsexed %>% filter(!is.na(sex))

all <- bind_rows(all, unsexed)

all %>% 
  as.data.frame() %>% 
  write_delim("complete_ARCHS4_age-sex_labels.tsv",
              delim = "\t", col_names = T)