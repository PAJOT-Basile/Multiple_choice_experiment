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

# Import the random SNPs that were selected
random_snps_selected <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Snps_selected_randomly.tsv",
                                   header = TRUE, sep = "\t")

# Metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

# Get the sizes of the chromosomes and compute where the next chromosome will start on the manhattan plot
chromosome_sizes <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Data/chromosome_sizes.tsv",
                               sep = "\t", header = TRUE)

# Import the sumstats of the thinned vcf (locus, chromosome, position, ...)
summary_stats <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Data/sumstats_thinned_50k.tsv",
                            sep = "\t", header = TRUE)

##############################
########### Run PCA ##########
##############################
pca <- scaleGen(data, NA.method="mean",scale=F,center=T) %>% 
  dudi.pca(scale=T, nf = 5,scannf = F)

# Extract the percentage of explained variance of interesting axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax4 <- ((pca$eig[4] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax5 <- ((pca$eig[5] / sum(pca$eig)) * 100) %>% round(digits = 2)

# Plot the PCA
pca$li %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(metadata, by = "ID_DNA_RAD") %>% 
  filter((Axis1 < -25 & Axis2 > 50) | Mother_ID == "M22-1") %>% 
  # mutate(toto = grepl("5300", ID_DNA_RAD)) %>% 
  ggplot(aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = Mother_ID)) +
  geom_text(aes(label = Label))


# Look at the contribution of each SNP along the genome to the axes of the PCA
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



