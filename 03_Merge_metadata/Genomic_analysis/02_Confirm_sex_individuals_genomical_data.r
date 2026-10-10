
##################################
#        WARNING
##################################
# This script was never used fully, so beware what you use it for.
























##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet", "patchwork", "ggforce", "ggnewscale", "ggh4x")
pacman::p_load(char = libraries, character.only = TRUE)

########################################
########### Useful functions ###########
########################################
source("98_Transversal_scripts/Transversal_functions.r")

########################################
########### Useful variables ###########
########################################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

colours_species <- c("praehirsuta" = "navy",
                     "forsmani" = "#3A9AB2",
                     "hybrid" = "red")

set.seed(123456789)

##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("Associated_data/raw_metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the extraction and library preparation duplicates
  filter(!grepl("_0", ID_DNA_RAD),
         # Remove individuals that were not fixed
         !is.na(ID_DNA_RAD)) %>% 
  rowwise() %>% 
  # Simplify the ID_DNA_RAD column
  mutate(ID_DNA_RAD_new = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all("a|b"),
         Family_level = factor(Family_level, levels = c("parents", "offspring")),
         Phenotype = factor(Phenotype, levels = c("forsmani", "hybrid", "praehirsuta")))

# Vcf
data <- read.vcfR("Associated_data/VCF_File.vcf") %>% 
  vcfR2genind()

# Some individuals are replicated (genotyped several times), so when we can, we 
# want to choose indivdiuals that have the less missing data.
## Find the replicates that have less missing data
replicate_to_keep <- (data[grepl("r|a|b", indNames(data))]@tab %>%
                        is.na() %>% rowSums() / (ncol(data@tab))) %>%
  as.data.frame() %>%
  rownames_to_column("ID_DNA_RAD") %>%
  rename(Missing_data = ".") %>% 
  mutate(ID_DNA_RAD_new = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all("a|b")) %>% 
  group_by(ID_DNA_RAD_new) %>% 
  filter(Missing_data == min(Missing_data)) %>% 
  pull(ID_DNA_RAD)

indivs_to_keep <- metadata %>% 
  filter(!grepl("r|a|b", ID_DNA_RAD)) %>% 
  pull(ID_DNA_RAD) %>% 
  c(replicate_to_keep)

# Subset the genind object
data <- data[indNames(data) %in% indivs_to_keep]

# Rename the individuals in the genind object to have simpler names
indNames(data) <- indNames(data) %>% 
  as.data.frame() %>% 
  rename(ID_DNA_RAD = ".") %>% 
  left_join(metadata, by = "ID_DNA_RAD") %>%
  pull(ID_DNA_RAD_new)

# Remove the duplicated columns in the metadata
metadata <- metadata %>% 
  select(-ID_DNA_RAD) %>% 
  rename(ID_DNA_RAD = ID_DNA_RAD_new) %>% 
  relocate(ID_DNA_RAD) %>% 
  unique()

summary_stats <- read.table("Associated_data/sumstats_thinned_50k.tsv",
                            header = TRUE, sep = "\t")

###################################################################
######### Find SNPs associated with sex in the two species ########
###################################################################
# To be sure that we get SNPs associated with sex, we need to keep data only for
# individuals for which we are sure about the sex. To do this, we keep only parents
# of the species. In addition, we need to use the output of the genomic PCA to 
# identify females from the different species
load("Associated_data/Rdata/Correspondence_name_species.rdata")

# As we know nothing about sexual determinism system in forsmani, individuals of the
# two species will be treated separately.


Hobs_species <- lapply(c("forsmani", "praehirsuta"), function(species, metadata, genetic_data, Species_names_parents){
  # First sample the genomic data to keep only parents of the praehirsuta species
  parents <- metadata %>%
    left_join(Species_names_parents, by = "ID_DNA_RAD") %>%
    filter(Family_level == "parents",
           Species == species) %>%
    pull(ID_DNA_RAD)

  data_sp <- genetic_data[indNames(genetic_data) %in% parents]

  data_sp@pop <- indNames(data_sp) %>%
    as.data.frame() %>%
    rename(ID_DNA_RAD = ".") %>%
    left_join(metadata %>%
                select(ID_DNA_RAD, Sex),
              by = "ID_DNA_RAD") %>%
    mutate(Sex = factor(Sex, levels = c("F", "M"))) %>%
    pull(Sex)


  Hobs_males <- compute_Hobs(data_sp[data_sp@pop == "M"])
  Hobs_females <- compute_Hobs(data_sp[data_sp@pop == "F"])

  Hobs_females %>%
    cbind(Hobs_males) %>%
    as.data.frame() %>%
    rownames_to_column("Position") %>%
    rename(Hobs_females = ".") %>%
    mutate(Locus_name = str_split_fixed(Position, "\\+|-", 2)[, 1] %>%
             str_remove_all(., ":$")) %>%
    left_join(summary_stats, by = "Locus_name") %>%
    mutate(Species = species) %>%
    return()

}, metadata, data, Species_names_parents) %>%
  bind_rows()

Hobs_plot <- Hobs_species %>% 
  select(-Position) %>% 
  unique() %>% 
  ggplot(aes(x = Hobs_males, y = Hobs_females)) + 
  geom_point() +
  facet_wrap(vars(Species), ncol = 1) +
  my_theme


Hobs_plot %>% 
  ggsave("../Figures/Suppfigs/Hobs_per_sex_species.png", .,
         width = 750, height = 550, units = "px",
         scale = 5)


# Individuals of the Jaera albifrons complex have a ZW sex determinism, so females
# are heterogametic and males are homogametic. To find SNPs associated with sex,
# we need to select SNPs that hare heterozygote in females and homozygotes in males
Sex_deterministic_SNPs <- Hobs_species %>% 
  filter(Hobs_males <= 0.05, Hobs_females >= 0.95) %>%
  select(Species, Position)

# Then, we run a PCA on these SNPs for the two species separately
pca_sex_determ_snps <- lapply(c("forsmani", "praehirsuta"), function(species, sex_determ_snps, metadata, genetic_data, Species_names_parents){
  parents <- metadata %>% 
    left_join(Species_names_parents, by = "ID_DNA_RAD") %>% 
    filter(Species == species) %>% 
    pull(ID_DNA_RAD)
  
  list_pos <- sex_determ_snps %>% 
    filter(Species == species) %>% 
    pull(Position) %>% 
    str_remove_all("\\.0|\\.1") %>%
    unique()
  
  data_filt <- genetic_data[loc = list_pos][indNames(genetic_data) %in% parents]
  print(data_filt)
  
  pca <- scaleGen(data_filt, NA.method = "mean", scale = F, center = T) %>%
    dudi.pca(scale = TRUE, nf = 5, scannf = FALSE)
  
  pca$li %>% 
    rownames_to_column("ID_DNA_RAD") %>% 
    left_join(metadata, by = "ID_DNA_RAD") %>% 
    mutate(Sex_pheno = case_when(Sex == "F" ~ "Females",
                                 Phenotype == "forsmani" ~ "Males forsmani",
                                 Phenotype == "praehirsuta" ~ "Males praehirsuta") %>%
             factor(levels = c("Males forsmani", "Males praehirsuta", "Females"))) %>% 
    return()
}, Sex_deterministic_SNPs, metadata, data, Species_names_parents) %>% 
  bind_rows()

pca_sex_determ_snps %>%
  left_join(Species_names_parents, by = "ID_DNA_RAD") %>% 
  ggplot(aes(x = Axis1, y = Axis2, colour = Sex_pheno, shape = Sex_pheno)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_colour_manual(name = "Sex and\nPhenotype",
                      values = c("Males forsmani" = "#3A9AB2",
                                 "Males praehirsuta" = "navy",
                                 "Females" = "grey50")) +
  scale_shape_manual(name = "Sex and\nPhenotype",
                     values = c("Males forsmani" = 17,
                                "Males praehirsuta" = 17,
                                "Females" = 15)) +
  facet_wrap(vars(Species)) +
  labs(x = "Axis 1",
       y = "Axis 2") +
  my_theme





# TODO: run the pca on LG3/LG7 on forsmani and on LG1 on praehirsuta
# TODO: tester pour voir si c'est bien conservé chez les descendants en utilisant l'espèce des parents comme indicateur

pca_sex_chroms <- lapply(c("forsmani", "praehirsuta"), function(species, list_chromosomes, metadata, genetic_data, summary_stats, Species_names_parents){
  lapply(list_chromosomes, function(chrom, species, metadata, genetic_data, summary_stats, Species_names_parents){
    parents <- metadata %>% 
      left_join(Species_names_parents, by = "ID_DNA_RAD") %>% 
      filter(
        # Family_level == "parents",
        Species == species) %>% 
      pull(ID_DNA_RAD)
    
    list_pos <- summary_stats %>% 
      filter(Chrom == chrom) %>% 
      mutate(Locus_name = paste0(Locus_name, ":+")) %>% 
      rbind(., 
            mutate(., Locus_name = str_replace_all(Locus_name, "\\+", "-"))) %>% 
      pull(Locus_name)
    
    data_filt <- genetic_data[loc = list_pos][indNames(genetic_data) %in% parents]
    
    pca <- scaleGen(data_filt, NA.method = "mean", scale = F, center = T) %>%
      dudi.pca(scale = TRUE, nf = 5,scannf = FALSE)
    
    pca$li %>% 
      rownames_to_column("ID_DNA_RAD") %>% 
      left_join(metadata, by = "ID_DNA_RAD") %>% 
      mutate(Sex_pheno = case_when(Sex == "F" ~ "Females",
                                   Phenotype == "forsmani" ~ "Males forsmani",
                                   Phenotype == "praehirsuta" ~ "Males praehirsuta") %>%
               factor(levels = c("Males forsmani", "Males praehirsuta", "Females"))) %>% 
      mutate(Chrom = chrom) %>% 
      return()
    
  }, species, metadata, genetic_data, summary_stats, Species_names_parents) %>% 
    bind_rows()
}, c("LG1", "LG3", "LG7"), metadata, data, summary_stats, Species_names_parents) %>% 
  bind_rows()



pca_sex_chroms %>% 
  filter(Chrom == "LG7") %>% 
  left_join(Species_names_parents, by = "ID_DNA_RAD") %>% 
  ggplot(aes(x = Axis1, y = Axis2, colour = Sex_pheno, shape = Sex_pheno)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_point(data = . %>% 
               filter(ID_DNA_RAD == "E4776"),
             size = 3, alpha = 0.7, colour = "red") +
  scale_colour_manual(name = "Sex and\nPhenotype",
                      values = c("Males forsmani" = "#3A9AB2",
                                 "Males praehirsuta" = "navy",
                                 "Females" = "grey50")) +
  scale_shape_manual(name = "Sex and\nPhenotype",
                     values = c("Males forsmani" = 17,
                                "Males praehirsuta" = 17,
                                "Females" = 15)) +
  facet_grid2(Species ~ Chrom, scales = "free") +
  labs(x = "Axis 1",
       y = "Axis 2") +
  my_theme
