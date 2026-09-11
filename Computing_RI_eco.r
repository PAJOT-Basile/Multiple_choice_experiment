##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = c("tidyverse"), character.only = TRUE)

##############################
###### Useful variables ######
##############################
set.seed(123456789)

nb_boots <- 10000

##############################
###### Useful functions ######
##############################
"%not_in%" <- function(x, y){!(x %in% y)}

compute_RI_in_bootstrap <- function(data_table, expected_substrates=c("praehirsuta_ascophyllum", "forsmani_pebble"), direction=NULL, resample=FALSE){
  if(resample) data_table <- resample_table(data_table)
  U <- dim(data_table[which(data_table$Phenotype_substrate %in% expected_substrates), ])[1]
  if (!is.null(direction)){
    S <- dim(data_table[which((data_table$Phenotype_substrate %not_in% expected_substrates) & (data_table$Phenotype == direction)), ])[1]
  }else{
    S <- dim(data_table[which(data_table$Phenotype_substrate %not_in% expected_substrates), ])[1]
  }
  
  return(1 - (S/(S+U)))
}

resample_table <- function(data_table, frac = 1){
  return(sample_n(data_table, size = round(frac * nrow(data_table)), replace = TRUE))
}

bootstrap_RI_eco <- function(data_table, direction=NULL, nb_bootstraps = nb_boots){
  boots_results <- replicate(nb_bootstraps,
                             expr = {compute_RI_in_bootstrap(data_table, direction=direction, resample=TRUE)},
                             simplify = "vector")
  return(boots_results)
}

##############################
######### Import data ########
##############################
sampling_data <- read.table("../../../Chapter_1_Sexual_selection/Code/Associated_data/Computing_RI_eco.tsv", sep = "\t", header = T) %>% 
  mutate(Phenotype_substrate = paste(Phenotype, Substrate, sep = "_"))


#############################################
######### Compute bootstraps on data ########
#############################################
# Praehirsuta
compute_RI_in_bootstrap(sampling_data, direction="praehirsuta")
bootstrap_results <- bootstrap_RI_eco(sampling_data, direction="praehirsuta")

# 95% Confidence interval
quantile(bootstrap_results, c(0.025, 0.975)) %>% 
  round(digits = 2)

# Forsmani
compute_RI_in_bootstrap(sampling_data, direction="forsmani")
bootstrap_results <- bootstrap_RI_eco(sampling_data, direction="forsmani")

# 95% Confidence interval
quantile(bootstrap_results, c(0.025, 0.975)) %>% 
  round(digits = 2)

# Both
compute_RI_in_bootstrap(sampling_data)
bootstrap_results <- bootstrap_RI_eco(sampling_data)

# 95% Confidence interval
quantile(bootstrap_results, c(0.025, 0.975)) %>% 
  round(digits = 2)

