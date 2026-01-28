##############################
########## Libraries #########
##############################
suppressMessages(if (!require("pacman")) install.packages("pacman"))
libraries <- c("tidyverse", "vcfR", "adegenet", "yaml")
suppressMessages(pacman::p_load(char = libraries, character.only = TRUE))
library(tidyverse); library(vcfR); library(adegenet); library(yaml); library(parallel)
cat("Making table marker\n")
##############################
########### Useful ###########
##############################
# Function not in
"%!in%" <- function(x, y){ return(!(x %in% y))}

# Compute_Error_Rate <- function(x, position, df){
#   simp_df <- df %>% 
#     filter(Sample_ID == x) %>% 
#     select(position)
#   
#   if (nrow(simp_df) < 2){
#     error_rate <- 0
#   }else{
#     differences_between_replicates <- sum(simp_df[1, ] != simp_df[2, ], na.rm = TRUE)
#     
#     error_rate <- differences_between_replicates 
#   }
#   return(error_rate)
# }
# 
# Get_False_Rate <- function(position, df){
#   list_samples <- df %>% 
#     pull(Sample_ID) %>% 
#     unique()
#   
#   error_rate <- sapply(list_samples, FUN = Compute_Error_Rate, position, df)
#   
#   return(sum(error_rate) / length(list_samples))
#   
# 
# }
# 
# Get_Error_Rate_Per_SNP <- function(df){
#   list_positions <- df %>% 
#     select(-Sample_ID) %>% 
#     names
#   
#   error_rate <- sapply(list_positions, FUN = Get_False_Rate, df)
#   return(error_rate)
# }


########################################


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
######### Import data ########
##############################

# Thinned vcf
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_and_without_missing_indivs.vcf",
                  verbose = FALSE) %>% 
  vcfR2genind()

# Import metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

# The order in which to arrange the SNPs
order_loci <- data@tab %>% 
  colnames() %>% 
  as_tibble() %>% 
  separate_wider_delim(value, names = c("Locus", "Col", "Allele"), ":")

# Get the sizes of the chromosomes and compute where the next chromosome will start on the manhattan plot
chromosome_sizes <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Data/chromosome_sizes.tsv",
                               sep = "\t", header = TRUE)

# Import the sumstats of the thinned vcf (locus, chromosome, position, ...)
summary_stats <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Data/sumstats_thinned_50k.tsv",
                            sep = "\t", header = TRUE)

##############################
######## Marker names ########
##############################
# Get the mean error rate using duplicates
SNPs_analysis <- data@tab %>% 
  as_tibble() %>% 
  mutate(sample = rownames(data@tab)) %>% 
  filter(grepl("r", sample)) %>% 
  select(contains(".0"), sample) %>% 
  mutate(Sample_ID = str_split_fixed(sample, "r", 2)[, 1]) %>% 
  select(-sample)

Error_rate <- Get_Error_Rate_Per_SNP(SNPs_analysis) %>%
  as_tibble() %>%
  rename(Error = value) %>%
  mutate(Position = names(Error_rate)) %>%
  relocate(Position)

Error_rate %>% 
  separate_wider_delim(Position, names = c("Locus", "Col", "Allele"), ":") %>% 
  filter(grepl("0", Allele)) %>% 
  mutate(Locus = Locus %>% as.integer,
         Locus_name = paste(Locus, Col, sep = ":"),
         Col = Col %>% as.integer) %>% 
  left_join(summary_stats, by = c("Locus", "Col", "Locus_name")) %>% 
  ggplot(aes(x = BP_cumul, y = Error)) +
  geom_point()

Error_rate %>%
  write.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Data/Error_rate_per_SNP.tsv",
              col.names = TRUE, row.names = FALSE, quote = FALSE)

data@tab %>%
  colnames() %>% 
  as_tibble() %>% 
  filter(grepl("\\.0", value)) %>% 
  mutate(value = str_remove_all(value, "\\+\\.0|\\-\\.0"),
         marker_type = 0,
         dropout = "0.0000",
         error = Error_rate) %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(toto = "",
         tata = c("!Marker names", "!Marker types, 0/1 = codominant/dominant", "!Allelic dropout rate", "!false allele rate")) %>% 
  relocate(toto) %>% 
  write.table(paste0(config_file$Outfolder, "colony.dat"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)