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
set.seed(1234)

all_snps <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/sumstats.tsv",
                   header = TRUE, sep = "\t")


random_snps_selected <- all_snps %>% 
  select(Locus.ID, Chr, BP, Col) %>% 
  mutate(Locus_name = paste(Locus.ID, Col, sep = ":")) %>% 
  group_by(Chr, Locus.ID) %>% 
  sample_n(size = 1)

# 
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


##############################
######### Import data ########
##############################
# random_snps_selected <- random_snps_selected %>% 
#   rename(Chromosome = Chr,
#          Locus = Locus.ID) %>% 
#   mutate(Chromosome = str_remove_all(Chromosome, "\\.1"))

random_snps_selected <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Snps_selected_randomly.tsv",
                                   header = TRUE, sep = "\t")
# Thinned vcf
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned.recode.vcf") %>% 
  vcfR2genind()

data <- data[loc = (locNames(data) %>% str_remove_all(":\\+|:\\-")) %in% random_snps_selected$Locus_name]

# Metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

chromosome_sizes <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/catalog.chrs.tsv",
                               sep = "\t", header = TRUE, comment = "") %>% 
  rename(Chromosome = X..Chrom) %>% 
  arrange(desc(Length)) %>%
  mutate(Chromosome = str_remove_all(Chromosome, "\\.1"),
         Cumul_Length = lag(cumsum(Length), default = 0))

summary_stats <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/populations.sumstats.tsv",
                            sep = "\t", header = TRUE,
                            skip = 1,
                            comment = "") %>% 
  rename(Locus = X..Locus.ID,
         Chromosome = Chr) %>% 
  select(Locus, Chromosome, BP, Col) %>% 
  unique %>% 
  mutate(Chromosome = str_remove_all(Chromosome, "\\.1")) %>% 
  left_join(chromosome_sizes %>% 
              select(Chromosome, Cumul_Length), by = "Chromosome") %>% 
  mutate(BP_cumul = BP + Cumul_Length) %>% 
  select(-Cumul_Length) %>% 
  right_join(random_snps_selected, by = c("Locus", "Chromosome", "BP", "Col"))


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
  Plot_Missing_Data_with_Threshold(0.5, dimension = "Indivs")

# See which individuals have the most missing data according to the threshold
missing_values_per_indiv %>% 
  filter(Proportion_missing_per_indiv >= 0.21) %>% 
  arrange(desc(Proportion_missing_per_indiv)) %>%
  left_join(metadata, by = join_by("Indiv" == "ID_DNA_RAD")) %>% 
  head(200) %>% 
  select(Indiv, Proportion_missing_per_indiv, Family_level, Sex)


##############################
#### Remove individuals with more than  50% of missing data####
##############################
indiv_less_than_50missing_data <- missing_values_per_indiv %>% 
  filter(Proportion_missing_per_indiv < 0.5) %>% 
  pull(Indiv)

# metadata %>% 
#   filter(ID_DNA_RAD %in% indiv_less_than_50missing_data) %>% 
#   pull(Label)

data_filt <- data[indiv_less_than_50missing_data]


pca <- scaleGen(data_filt, NA.method="mean",scale=F,center=T) %>% 
  dudi.pca(scale=T, nf = 5,scannf = F)

# Extract the percentage of explained variance of interesting axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax4 <- ((pca$eig[4] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax5 <- ((pca$eig[5] / sum(pca$eig)) * 100) %>% round(digits = 2)

pca$li %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(metadata, by = "ID_DNA_RAD") %>% 
  mutate(toto = grepl("5300", ID_DNA_RAD)) %>% 
  ggplot() +
  geom_point(aes(x = Axis1, y = Axis2, color = Phenotype))



pca$co %>% 
  rownames_to_column("Position") %>% 
  separate_wider_delim(Position, names = c("Locus", "Col", "Allele"), ":") %>% 
  filter(grepl("0", Allele)) %>% 
  mutate(Locus = Locus %>% as.integer,
         Locus_name = paste(Locus, Col, sep = ":"),
         Col = Col %>% as.integer) %>% 
  left_join(summary_stats, by = c("Locus", "Col", "Locus_name")) %>% 
  ggplot(aes(x = BP_cumul, y = abs(Comp1))) +
  geom_point()


