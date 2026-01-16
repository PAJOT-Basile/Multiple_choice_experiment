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

# Function not in
"%!in%" <- function(x, y){ return(!(x %in% y))}

##############################
######### Import data ########
##############################
# Thinned vcf
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_and_without_missing_indivs.vcf") %>% 
  vcfR2genind()

random_snps_selected <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Snps_selected_randomly.tsv",
                                   header = TRUE, sep = "\t")

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



