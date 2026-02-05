##############################
########## Libraries #########
##############################
suppressMessages(if (!require("pacman")) install.packages("pacman"))
libraries <- c("tidyverse", "vcfR", "adegenet", "yaml", "parallel")
suppressMessages(pacman::p_load(char = libraries, character.only = TRUE))

cat("Making table marker\n")
##############################
########### Useful ###########
##############################
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
######### Import data ########
##############################
# Import the parameter file
config_file <- read_yaml("Data/Input.yaml")
# Thinned vcf
data <- read.vcfR(config_file$Input_file, verbose = FALSE) %>% 
  vcfR2genind()

# Import metadata
metadata <- read.table(config_file$Metadata,
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

##############################
######## Marker names ########
##############################
# Get the mean error rate using duplicates
suppressWarnings(SNPs_analysis <- data@tab %>% 
  as_tibble() %>% 
  mutate(sample = rownames(data@tab)) %>% 
  filter(grepl("r", sample)) %>% 
  select(contains(".0"), sample) %>% 
  mutate(Sample_ID = str_split_fixed(sample, "r", 2)[, 1]) %>% 
  select(-sample))

Error_rate <- Get_Error_Rate_Per_SNP(SNPs_analysis)
suppressWarnings( Error_rate %>%
  as_tibble() %>%
  rename(Error = value) %>%
  mutate(Position = names(Error_rate),
         Error = round(Error, digits = 4)) %>%
  mutate(Position = str_remove_all(Position, "\\+\\.0|\\-\\.0"),
         marker_type = 0, dropout_rate = 0) %>%
  relocate(Position, marker_type, dropout_rate, Error) %>%
  t() %>%
  as_tibble() %>%
  mutate(toto = "",
  tata = c("!Marker names", "!Marker types, 0/1 = codominant/dominant", "!Allelic dropout rate", "!false allele rate")) %>% 
  relocate(toto) %>% 
  write.table(paste0(config_file$Outfolder, "colony.dat"),
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE))