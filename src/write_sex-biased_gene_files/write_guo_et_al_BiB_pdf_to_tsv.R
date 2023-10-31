# turn pdf of results from https://academic.oup.com/bib/article/19/2/188/2731023#120301544
# into tsv
library(tidyverse)
library(pdftools)

# read in pdf, which is a list with one element per page
guo <- pdftools::pdf_text(pdf = "~/projects/age-sex-prediction/src/write_sex-biased_gene_files/guo_supp.pdf")

# string split on new line (makes nested lists)
guo <- lapply(guo, str_split, pattern = "\n")
# now each element in the list is a character vetor instead of list
guo <- lapply(guo, unlist)
# make one character vector
guo <- unlist(guo)
# remove title  from first page entry
guo <- guo[guo != "supplementary table S1. Sex-biased genes by tissues"]
# get rid of empty strings
guo <- guo[guo != ""]
# get rid of lines with this url from each page
guo <- guo[!str_detect(guo, "http://mc.manuscriptcentral.com/bib")]
# turn to tibble
guo <- as_tibble(guo)
# get rid of column names
guo <- guo[-1,]
# separate cols on white space
# use merge so it only uses first white space
# and keeps rest of items
guo <- guo %>% 
  separate(value, into = c("hgnc_symbol", "values"), sep = "\\s+", extra = "merge") %>% 
  separate(values, into = c("chromosome", "values"), sep = "\\s+", extra = "merge") %>% 
  separate(values, into = c("sex_bias", "tissue"), sep = "\\s+", extra = "merge")

guo %>% 
  write_delim("~/projects/age-sex-prediction/src/write_sex-biased_gene_files/guo_supp.tsv",
              delim = "\t", col_names = T)

