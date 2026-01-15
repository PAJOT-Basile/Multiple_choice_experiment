##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
##### Select random snps #####
##############################
set.seed(1234)

all_snps <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/sumstats.tsv",
                   header = TRUE, sep = "\t")


x <- all_snps %>% 
  select(Locus.ID, Chr, BP) %>% 
  group_by(Chr, Locus.ID) %>% 
  sample_n(size = 1)


x %>%
  ungroup %>%
  select(Chr, BP) %>%
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/whitelist.tsv",
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rm(list = ls())
##############################
##### Import thinned vcf #####
##############################

data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/thinned.vcf") %>% 
  vcfR2genind()
