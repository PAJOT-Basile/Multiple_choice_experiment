##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
########### Useful ###########
##############################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

# Function to plot the distribution of missing values and count how many SNPs are 
# removed when given a threshold
Plot_Missing_Data_with_Threshold <- function(df, threshold_missing_data = 0.0, dimension = "SNP"){
  # Compute the number of SNPs that will be removed
  removed <- nrow(df) - (df %>% 
    filter(Proportion_missing <= threshold_missing_data) %>% 
    nrow)
  
  # Look at the distribution of missing data
  return(df %>% 
           ggplot(aes(x = Proportion_missing)) +
           geom_histogram(bins = 200) +
           geom_vline(xintercept = threshold_missing_data, color = "red") +
           labs(subtitle = paste0(threshold_missing_data * 100, "% of NA: remove ", removed, " ", dimension, " (",
                                  round(removed / nrow(df), digits = 2) * 100, "%)."),
                y = paste0(str_remove_all(dimension, "s"), " count"),
                x = paste0("Proportion of missing data per ", str_remove_all(dimension, "s"))) +
           my_theme)
}


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


##############################
#### Analyse missing data ####
##############################
# Compute the proportion of missing values per SNP
missing_values_per_snp <- data@tab %>% 
  is.na %>% 
  colSums %>% 
  as_tibble %>% 
  rename(Nb_missing_values = value) %>% 
  mutate(SNP = colnames(data@tab),
         Proportion_missing_per_SNP = Nb_missing_values / nrow(data@tab)) %>% 
  relocate(SNP) %>% 
  separate_wider_delim(SNP, ":", names = c("Locus", "Position_on_Locus", "Allele")) %>% 
  filter(grepl("0", Allele)) %>% 
  unite(SNP, c(Locus, Position_on_Locus, Allele), sep = ":")

missing_values_per_snp %>% 
  rename(Proportion_missing = Proportion_missing_per_SNP) %>% 
  Plot_Missing_Data_with_Threshold(0.1)

# Do the same for the proportion of missing values per individual
missing_values_per_indiv <- data@tab %>% 
  is.na %>% 
  rowSums %>% 
  as_tibble %>% 
  rename(Nb_missing_values = value) %>% 
  mutate(Indiv = rownames(data@tab),
         Proportion_missing_per_indiv = Nb_missing_values / ncol(data@tab)) %>% 
  relocate(Indiv)

missing_values_per_indiv %>% 
  rename(Proportion_missing = Proportion_missing_per_indiv) %>% 
  Plot_Missing_Data_with_Threshold(0.25, dimension = "Indivs")
