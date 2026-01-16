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
Plot_Missing_Data_with_Threshold <- function(df, threshold_missing_data = 0.0, dimension = "SNP", ...){
  # Compute the number of SNPs that will be removed
  removed <- nrow(df) - (df %>% 
    filter(Proportion_missing <= threshold_missing_data) %>% 
    nrow)
  
  # Look at the distribution of missing data
  return(df %>% 
           ggplot(aes(x = Proportion_missing)) +
           geom_histogram(...) +
           geom_vline(xintercept = threshold_missing_data, color = "red") +
           labs(subtitle = paste0(threshold_missing_data * 100, "% of NA: remove ", removed, " ", dimension, " (",
                                  round(removed / nrow(df), digits = 2) * 100, "%)."),
                y = paste0(str_remove_all(dimension, "s"), " count"),
                x = paste0("Proportion of missing data per ", str_remove_all(dimension, "s"))) +
           my_theme)
}

# Function not in
"%!in%" <- function(x, y){ return(!(x %in% y))}

##############################
##### Select random snps #####
##############################
# Set the random seed to make the random selection of SNPs replicable
set.seed(1234)

# Get the list of all SNPs that were outputed from populations
all_snps <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/sumstats.tsv",
                   header = TRUE, sep = "\t")


# Select one random SNP per locus and per chromosome
random_snps_selected <- all_snps %>% 
  select(Locus.ID, Chr, BP, Col) %>% 
  mutate(Locus_name = paste(Locus.ID, Col, sep = ":")) %>% 
  group_by(Chr, Locus.ID) %>% 
  sample_n(size = 1)

# Uncomment to filter the vcf
random_snps_selected %>%
  ungroup %>%
  select(Chr, BP) %>%
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/whitelist.tsv",
                  sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
        args = paste0(" --vcf /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/populations.snps.vcf",
                      " --out /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned",
                      " --positions /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/whitelist.tsv",
                      " --recode"))


#####################
#### SNP removal ####
#####################
# Some locus are present several times on the same position in the genome. Thus,
# when running vcftools, the filtering process keeps the same SNP several times
# (from the different RADtags at the same position). In this part, we filter out
# the unwanted SNPs.

# Import the vcf file to check the list of loci
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned.recode.vcf") %>% 
  vcfR2genind()

# Look for loci that are not wanted in the vcf file
unwanted_snps <- locNames(data)[(str_remove_all(locNames(data), ":\\+|:\\-") %!in% random_snps_selected$Locus_name)]

# Make a file containing the unwanted SNPs
unwanted_snps %>% 
  as_tibble %>% 
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/unwanted_snps.txt",
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Call grep to filter out the unwanted SNPs
system2("grep", args = paste0("-v -f /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/unwanted_snps.txt",
                              " /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned.recode.vcf",
                              " > /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps.vcf"))

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

# No problem on the SNPs, let's look at the individuals

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
  Plot_Missing_Data_with_Threshold(0.5, dimension = "Indivs", bins = 50)

# We have individuals with a lot of missing values, so we need to filter out
# some individuals to run the PCA
indiv_less_than_50missing_data <- missing_values_per_indiv %>% 
  filter(Proportion_missing_per_indiv < 0.5) %>% 
  select(Indiv)

indiv_less_than_50missing_data %>% 
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Indivs_to_keep.txt",
              sep = "t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Remove individuals that have too much missing data
system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
        args = paste0("--vcf /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps.vcf",
                      " --keep /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Indivs_to_keep.txt",
                      " --recode --stdout",
                      " > /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_and_without_missing_indivs.vcf"))