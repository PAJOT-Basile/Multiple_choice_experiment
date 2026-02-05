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


Get_Error_Rate_One_Position <- function(position, df){
  Compute_Error_Pair_Replicates <- function(pair, simp_df){
    genotypes <- simp_df[
      which(simp_df["Sample_ID"] == pair),
      setdiff(names(simp_df), "Sample_ID")
    ]
    if (length(genotypes) < 2) rate_to_return <- 0
    else{
      nb_differences <- sum(genotypes[1] != genotypes[2], na.rm = TRUE)
      nb_NA <- sum(is.na(genotypes))
      
      rate_to_return <- max(nb_differences - nb_NA, 0)
    }
    return(rate_to_return)
  }
  simp_df <- df[c("Sample_ID", position)]
  
  list_pairs <- unique(df$Sample_ID) 
  nb_pairs <- length(list_pairs)
  
  error_pair <- sapply(list_pairs, Compute_Error_Pair_Replicates, simp_df)
  return(sum(error_pair) / nb_pairs)
}

Get_Error_Rate_Per_SNP <- function(df){
  list_positions <- df %>% 
    select(-Sample_ID) %>% 
    names()
  # prepare parallel computations
  cl <- makeCluster(20)
  error_rate <- parSapply(cl = cl, X = list_positions, FUN = Get_Error_Rate_One_Position, df)
  stopCluster(cl)
  return(error_rate)
}

##############################
##### Select random snps #####
##############################
# Set the random seed to make the random selection of SNPs replicable
# set.seed(1234)
# 
# # Get the list of all SNPs that were outputed from populations
# all_snps <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/sumstats.tsv",
#                    header = TRUE, sep = "\t")
# 
# 
# # Select one random SNP per locus and per chromosome
# random_snps_selected <- all_snps %>% 
#   select(Locus.ID, Chr, BP, Col) %>% 
#   mutate(Locus_name = paste(Locus.ID, Col, sep = ":")) %>% 
#   group_by(Chr, Locus.ID) %>% 
#   sample_n(size = 1)
# 
# # Uncomment to filter the vcf
# random_snps_selected %>%
#   ungroup %>%
#   select(Chr, BP) %>%
#   write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/whitelist.tsv",
#                   sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
#         args = paste0(" --vcf /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/populations.snps.vcf",
#                       " --out /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned",
#                       " --positions /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/whitelist.tsv",
#                       " --recode"))
# 
# 
# #####################
# #### SNP removal ####
# #####################
# # Some locus are present several times on the same position in the genome. Thus,
# # when running vcftools, the filtering process keeps the same SNP several times
# # (from the different RADtags at the same position). In this part, we filter out
# # the unwanted SNPs.
# 
# # Import the vcf file to check the list of loci
# data <- read.vcfR("../../Manips/Manip_sexualselection1/Data/Genomic/thinned.vcf") %>% 
#   vcfR2genind()
# 
# # Look for loci that are not wanted in the vcf file
# unwanted_snps <- locNames(data)[(str_remove_all(locNames(data), ":\\+|:\\-") %!in% random_snps_selected$Locus_name)]
# 
# # Make a file containing the unwanted SNPs
# unwanted_snps %>% 
#   as_tibble %>% 
#   write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/unwanted_snps.txt",
#               sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# 
# # Call grep to filter out the unwanted SNPs
# system2("grep", args = paste0("-v -f /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/unwanted_snps.txt",
#                               " /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned.recode.vcf",
#                               " > /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps.vcf"))

##############################
######### Import data ########
##############################

# Thinned vcf
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps.vcf",
                  verbose = FALSE) %>% 
  vcfR2genind()

# Import metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

##############################
######## Marker names ########
##############################
# Get the mean error rate using duplicates
# SNPs_analysis <- data@tab %>% 
#   as_tibble() %>% 
#   mutate(sample = rownames(data@tab)) %>% 
#   filter(grepl("r", sample)) %>% 
#   select(contains(".0"), sample) %>% 
#   mutate(Sample_ID = str_split_fixed(sample, "r", 2)[, 1]) %>% 
#   select(-sample)
# 
# Error_rate <- Get_Error_Rate_Per_SNP(SNPs_analysis) %>%
#   as_tibble() %>%
#   rename(Error = value) %>%
#   mutate(Position = names(Error_rate)) %>%
#   relocate(Position)

Error_rate <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Error_rate_per_SNP.tsv",
                         sep = " ", header = TRUE)

# Filter to remove the SNPs with errors
Error_rate %>% 
  filter(Error > 0) %>% 
  select(Position) %>% 
  mutate(Position = str_split_fixed(Position, "\\.0|\\.1", 2)[, 1]) %>% 
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Positions_high_error.tsv",
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

system2("grep", args = paste0("-v -f /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Positions_high_error.tsv",
                              " /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps.vcf",
                              " > /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_no_error.vcf"))


# Import data with no error
data_no_err <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_no_error.vcf",
                         verbose = FALSE) %>% 
  vcfR2genind()
##############################
#### Analyse missing data ####
##############################
# Compute the proportion of missing values per SNP
missing_values_per_snp <- data_no_err@tab %>% 
  is.na %>% 
  colSums %>% 
  as_tibble %>% 
  rename(Nb_missing_values = value) %>% 
  mutate(SNP = colnames(data_no_err@tab),
         Proportion_missing_per_SNP = Nb_missing_values / nrow(data_no_err@tab)) %>% 
  relocate(SNP) %>% 
  separate_wider_delim(SNP, ":", names = c("Locus", "Position_on_Locus", "Allele")) %>% 
  filter(grepl("0", Allele)) %>% 
  unite(SNP, c(Locus, Position_on_Locus, Allele), sep = ":")

missing_values_per_snp %>% 
  rename(Proportion_missing = Proportion_missing_per_SNP) %>% 
  Plot_Missing_Data_with_Threshold(0.2, bins = 300)


# No problem on the SNPs, let's look at the individuals

# Do the same for the proportion of missing values per individual
missing_values_per_indiv <- data_no_err@tab %>% 
  is.na %>% 
  rowSums %>% 
  as_tibble %>% 
  rename(Nb_missing_values = value) %>% 
  mutate(Indiv = rownames(data_no_err@tab),
         Proportion_missing_per_indiv = Nb_missing_values / ncol(data@tab)) %>% 
  relocate(Indiv)

missing_values_per_indiv %>% 
  rename(Proportion_missing = Proportion_missing_per_indiv) %>% 
  Plot_Missing_Data_with_Threshold(0.75, dimension = "Indivs", bins = 50)

# We have individuals with a lot of missing values, so we need to filter them out
# In particular, the two extreme individuals were fixed dead and did not map particularly well, 
# so we will take them out here
# Make a list of individuals to keep
indiv_keep <- missing_values_per_indiv %>% 
  filter(Proportion_missing_per_indiv < 0.75) %>% 
  select(Indiv)

indiv_keep %>% 
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Indivs_to_keep_missing_data.txt",
              sep = "t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Remove individuals that have too much missing data
system2("/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools",
        args = paste0("--vcf /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_no_error.vcf",
                      " --keep /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Indivs_to_keep_missing_data.txt",
                      " --recode --stdout",
                      " > /shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_no_error_and_good_indivs.vcf"))
