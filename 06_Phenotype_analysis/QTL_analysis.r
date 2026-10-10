# Here, we want to try to map the standardised residuals of the number of spines
# and the number of curved setae. To do this, the steps are:
#     1. Import data
#     2. Separate sexes
#     3. Make a linear regression and extract the standardised residuals
#     4. Prepare all the data to run the QTL analysis (phentype matrix, 
#         relatedness matrix between individuals, genomic data)
#     5. Run the QTL analysis


########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "adegenet", "vcfR", "statgenGWAS", "ggforce", "ggh4x", "ggpubr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
###### Useful variables ######
##############################
set.seed(123456789)
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20))

theme_grid <- my_theme +
  theme(plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.margin = margin(t = 15, r = 15, b = 15, l = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15),
        text = element_text(size = 12))

vcf_file <- "/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/VCF_File.vcf"
vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"

group_names <- c("Male_forsmani", "Male_praehirsuta", "Female_forsmani", "Female_praehirsuta")

##############################
###### Useful functions ######
##############################
"%not_in%" <- function(x, y){return(!(x %in% y))}

Get_chromosome_summary_information <- function(chromosome_recap){
  #' This function makes a recap of the chromosome information
  #' (max position, center, ...).
  #' /!\ THE CHROMOSOMES HAVE TO BE IN THE CORRECT ORDER AND THEIR COLUMN HAS TO BE NAMED "Chromosome"
  if ("Chromosome" %not_in% names(chromosome_recap)) stop("Is your column containing the names of the chromosomes named 'Chromosome'?")
  if ("Position" %not_in% names(chromosome_recap)) stop("Is your column containing the positions along the chromosomes named 'Position'?")
  
  
  chromosome_recap %>% 
    group_by(Chromosome) %>% 
    summarize(max_position = max(Position, na.rm = TRUE)) %>% 
    mutate(Start_chromosome = cumsum(lag(max_position + 1, default = 0)),
           Center_chromosome = round(max_position) / 2 + Start_chromosome) %>% 
    return()
}


##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the dulicates
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "offspring") %>% 
  rowwise() %>% 
  mutate(ID_DNA_RAD1 = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all("a|b"),
         Size = str_replace_all(Size, ",", ".") %>% 
           as.numeric(),
         # Recompute the total number of spines and setae per pereopod to homogeneise the values for the parents and offspring
         P7.ep = case_when(
           all(is.na(c(P7.very.large.ep, P7.larg.ep, P7.med.ep, P7.small.ep))) ~ NA,
           TRUE ~ sum(P7.very.large.ep, P7.larg.ep, P7.med.ep, P7.small.ep, na.rm = TRUE)),
         P6.ep.dist = case_when(all(is.na(c(P6.ep.dist.very.large, P6.ep.dist.medium, P6.ep.dist.large, P6.ep.dist.small))) ~ NA,
                                TRUE ~ sum(P6.ep.dist.very.large, P6.ep.dist.medium, P6.ep.dist.large, P6.ep.dist.small, na.rm = TRUE)),
         P6.ep = case_when(
           all(is.na(c(P6.very.large.ep, P6.larg.ep, P6.med.ep, P6.small.ep))) ~ NA,
           TRUE ~ sum(P6.very.large.ep, P6.larg.ep, P6.med.ep, P6.small.ep, na.rm = TRUE)
         ),
         Family_level = factor(Family_level, levels = c("parents", "offspring")),
         Species = factor(Species, levels = c("forsmani", "praehirsuta")),
         # Remove the ">" that are in the data
         across(starts_with(paste0("P", 1:7)), ~ str_remove_all(., ">") %>% trimws() %>% as.numeric())) %>% 
  unique() %>% 
  mutate(Sum_curved_setaeP1P5 = case_when(is.na(P1.curv.setae) & is.na(P2.curv.setae) & is.na(P3.curv.setae) & is.na(P4.curv.setae) & is.na(P5.curv.setae) ~ NA,
                                          TRUE ~ sum(P1.curv.setae, P2.curv.setae, P3.curv.setae, P4.curv.setae, P5.curv.setae, na.rm = TRUE)),
         Sum_spinesP4P7 = case_when(is.na(P4.small.ep) & is.na(P5.ep) & is.na(P6.ep) & is.na(P7.ep) ~ NA,
                                    TRUE ~ sum(P4.small.ep, P5.ep, P6.ep, P7.ep, na.rm = TRUE))) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, Sum_curved_setaeP1P5, Sum_spinesP4P7), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()



data <- read.vcfR(vcf_file, verbose = FALSE) %>% 
  vcfR2genind()

# Replicates with less missing data
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
  select(ID_DNA_RAD) %>% 
  drop_na() %>% 
  pull(ID_DNA_RAD) %>% 
  c(replicate_to_keep)

# Subset the genind object
data <- data[indNames(data) %in% indivs_to_keep]
# Add the sex of the individuals to the genomic data to be able to separate them easily
data@other$Sex <- rownames(data@tab) %>% 
  as.data.frame() %>% 
  rename(ID_DNA_RAD = ".") %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Sex),
            by = "ID_DNA_RAD") %>% 
  pull(Sex) %>% 
  as.factor()
# Add the species of the individuals to separate them easily
data@other$Species <- rownames(data@tab) %>% 
  as.data.frame() %>% 
  rename(ID_DNA_RAD = ".") %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species),
            by = "ID_DNA_RAD") %>% 
  pull(Species) %>% 
  as.factor()

# Do the same in the metadata
metadata <- metadata %>% 
  filter(ID_DNA_RAD %in% indivs_to_keep)

# Define the four groups that we need to use
groups_individuals <- lapply(c("M", "F"), function(sex, df){
  lapply(c("forsmani", "praehirsuta"), function(species, sex, df){
    df %>% 
      filter(Sex == sex, Species == species) %>% 
      pull(ID_DNA_RAD)
  }, sex, df)
}, metadata) %>% 
  unlist(recursive = FALSE)
names(groups_individuals) <- group_names


# Import the sumstats of the thinned vcf (locus, chromosome, position, ...)
summary_stats <- read.table("Associated_data/sumstats_thinned_50k.tsv",
                            sep = "\t", header = TRUE)

######################################
############ Make the lms ############
######################################
# Here, a simple linear regression was tested as well as a quadratic regression
# between size and the number of curved setae/spines. The best model was selected
# using AIC.
metadata %>% 
  filter(Sex == "M", Species == "forsmani") %>% 
  ggplot(aes(x = Size, y = Sum_spinesP4P7)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x + I(x^2)) +
  geom_point()
# Males forsmani spines
lm_forsmani_spines <- metadata %>%
  filter(Sex == "M", Species == "forsmani") %>% 
  lm(Sum_spinesP4P7 ~ Size + I(Size^2), data = .)

# Males forsmani setae
lm_forsmani_setae <- metadata %>%
  filter(Sex == "M", Species == "forsmani") %>% 
  lm(Sum_curved_setaeP1P5  ~ Size, data = .)

# Males praehirsuta setae
lm_praehirsuta_setae <- metadata %>%
  filter(Sex == "M", Species == "praehirsuta") %>% 
  lm(Sum_curved_setaeP1P5 ~ Size, data = .)
  
# Males praehirsuta spines
lm_praehirsuta_spines <- metadata %>% 
  filter(Sex == "M", Species == "praehirsuta") %>% 
  lm(Sum_spinesP4P7 ~ Size, data = .)

lms <- list("forsmani_spines" = lm_forsmani_spines,
            "forsmani_setae" = lm_forsmani_setae,
            "praehirsuta_spines" = lm_praehirsuta_spines,
            "praehirsuta_setae" = lm_praehirsuta_setae)


######################################
######### Compute relatedness ########
######################################
# First, make lists of males and females from the two species to keep
# lapply(c("forsmani", "praehirsuta"), function(species, df, vcf_file){
#   print(species)
#   metadata %>%
#     filter(Species == species) %>%
#     select(ID_DNA_RAD) %>%
#     write.table(paste0("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Useful_vcftools/", species, "/Offspring.txt"),
#                 col.names = FALSE, quote = FALSE, row.names = FALSE)
# 
#   # Compute the relatedness between offspring
#   system2(vcftools, args = paste0("--vcf ", vcf_file,
#                                   " --relatedness",
#                                   " --keep ", "/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Useful_vcftools/", species, "/Offspring.txt",
#                                   " --out ", "Associated_data/relatedness_", species))
# }, metadata, vcf_file)

# Import the relatedness
relatedness <- lapply(c("forsmani", "praehirsuta"), function(species){
  table <- read.table(paste0("Associated_data/relatedness_", species, ".relatedness"), header = TRUE) %>%
    rename(Relatedness = RELATEDNESS_AJK) %>%
    mutate(Species = species)
    return(table)
}) %>%
  bind_rows()

# Plot the relatedness
relatedness %>%
  ggplot(aes(x = INDV1, y = INDV2, fill = Relatedness)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("navy", "dodgerblue", "skyblue", "white", "orange2", "red", "firebrick")) +
  facet_wrap(vars(Species), scales = "free")

# Scale the relatedness
scaled_relatedness <- lapply(c("Male", "Female"), function(sex, relatedness, metadata){
  lapply(c("forsmani", "praehirsuta"), function(species, sex, relatedness, metadata){
    tab <- relatedness %>% 
      filter(Species == species) %>% 
      left_join(metadata %>% 
                  select(ID_DNA_RAD, Sex, Size, Sum_curved_setaeP1P5, Sum_spinesP4P7),
                by = join_by("INDV1" == "ID_DNA_RAD")) %>% 
      left_join(metadata %>% 
                  select(ID_DNA_RAD, Sex, Size, Sum_curved_setaeP1P5, Sum_spinesP4P7),
                by = join_by("INDV2" == "ID_DNA_RAD")) %>% 
      filter(Sex.x == str_split(sex, "", simplify = TRUE)[, 1],
             Sex.y == str_split(sex, "", simplify = TRUE)[, 1],
             !is.na(Size.x), !is.na(Size.y),
             !is.na(Sum_curved_setaeP1P5.x), !is.na(Sum_curved_setaeP1P5.y),
             !is.na(Sum_spinesP4P7.x), !is.na(Sum_spinesP4P7.y)) %>% 
      select(-c(starts_with("Sex."), starts_with("Size"), starts_with("Sum"))) %>%
      select(-Species) %>%
      pivot_wider(names_from = INDV2, values_from = Relatedness) %>%
      column_to_rownames("INDV1")
    for (i in 1:nrow(tab)){
      for (j in 1:ncol(tab)){
        if (i == j){
          next
        }else if (tab[i, j] %>% is.na){
          tab[i, j] <- tab[j, i]
        }
      }
    }
    min_rel <- min(tab)
    max_rel <- max(tab)
    tab %>%
      mutate(across(everything(), .fns = function(x) (x + abs(min_rel)) / (max_rel - min_rel))) %>%
      as.matrix() %>%
      return()
  }, sex, relatedness, metadata)
}, relatedness, metadata) %>% 
  unlist(recursive = FALSE)
names(scaled_relatedness) <- group_names

###########################################
######### Prepare the genomic data ########
###########################################
# Scale the genome to get rid of missing data
X <- scaleGen(data, NA.method = "mean", scale = FALSE, center = TRUE)

map_chromosome <- colnames(X) %>% 
    as.data.frame() %>%
    rename(Position = ".") %>%
    mutate(pos = Position) %>%
    separate(pos, c("Locus", "Col", "Allele"), ":") %>%
    mutate(across(c(Col, Locus), ~ as.numeric(.))) %>%
    left_join(summary_stats, by = c("Locus", "Col")) %>%
    select(Position, Chrom, BP) %>%
    column_to_rownames("Position") %>% 
    rename(chr = Chrom, pos = BP)

geno <- X %>% 
    as.data.frame() %>% 
    rownames_to_column("ID_DNA_RAD") %>% 
    inner_join(metadata,
               by = "ID_DNA_RAD") %>% 
    select(ID_DNA_RAD, rownames(map_chromosome)) %>% 
    column_to_rownames("ID_DNA_RAD")

###########################################
######### Separate the phenotypes #########
###########################################
phenos <- lapply(c("Male", "Female"), function(sex, df){
  lapply(c("forsmani", "praehirsuta"), function(species, sex, df){
    df %>% 
      filter(Sex == str_split(sex, "", simplify = TRUE)[, 1],
             Species == species,
             !is.na(ID_DNA_RAD)) %>% 
      select(ID_DNA_RAD, Size, Sum_curved_setaeP1P5, Sum_spinesP4P7) %>% 
      select_if(~sum(!is.na(.)) > 0) %>% 
      mutate(genotype = ID_DNA_RAD) %>% 
      column_to_rownames("ID_DNA_RAD") %>% 
      relocate(genotype) %>% 
      drop_na()
  }, sex, df)
}, metadata) %>% 
  unlist(recursive = FALSE)
# phenos <- lapply(c("Male", "Female"), function(sex, df, lms){
#   lapply(c("forsmani", "praehirsuta"), function(species, df, lms, sex){
#     if (sex == "Male"){
#       lapply(c("setae", "spines", "size"), function(trait, lms, df, species){
#         if (trait != "size"){
#           to_ret <- lms[[paste0(species, "_", trait)]]$residuals %>% 
#             as.data.frame() %>% 
#             rename(!!sym(paste0("Sum_", trait)) := ".")
#         }else{
#           to_ret <-  metadata %>% 
#             filter(Sex == str_split(sex, "", simplify = TRUE)[, 1],
#                    Species == species,
#                    !is.na(Size)) %>% 
#             select(Size)
#         }
#       }, lms, df, species) %>% 
#         bind_cols() %>%
#         cbind(
#           metadata %>% 
#             filter(Sex == str_split(sex, "", simplify = TRUE)[, 1],
#                    Species == species,
#                    !is.na(Size)) %>% 
#             select(ID_DNA_RAD)) %>% 
#         as_tibble() %>% 
#         rename(genotype = ID_DNA_RAD) %>% 
#         mutate(toto = genotype) %>% 
#         column_to_rownames("toto") %>% 
#         relocate(genotype)
#     }else{
#       metadata %>% 
#         filter(Sex == str_split(sex, "", simplify = TRUE)[, 1],
#                Species == species,
#                !is.na(Size)) %>% 
#         select(ID_DNA_RAD, Size) %>% 
#         mutate(toto = ID_DNA_RAD) %>% 
#         rename(genotype = ID_DNA_RAD) %>% 
#         column_to_rownames("toto") %>% 
#         relocate(genotype)
#     }
#   }, df, lms, sex)
# }, metadata, lms) %>% 
#   unlist(recursive = FALSE)
names(phenos) <- group_names

###########################################
############## Run the GWAS ###############
###########################################
gwas_out <- lapply(group_names, function(group_name, geno, map_chromosome, phenos, scaled_relatedness, df){
  print(group_name)
  phenos_group <- phenos[[group_name]]
  gdata_object <- createGData(geno = geno,
              map = map_chromosome,
              pheno = phenos_group,
              kin = scaled_relatedness[[group_name]])

  runSingleTraitGwas(gdata_object, thrType = "bonf")
  
}, geno, map_chromosome, phenos, scaled_relatedness, metadata)
names(gwas_out) <- group_names


GWAS_out <- gwas_out$Male_forsmani$GWAResult$phenos_group %>% 
  filter(!is.na(pValue)) %>% 
  mutate(Sex = "Male", Species = "forsmani") %>% 
  rbind(gwas_out$Male_praehirsuta$GWAResult$phenos_group %>% 
          filter(!is.na(pValue)) %>% 
          mutate(Sex = "Male", Species = "praehirsuta")) %>% 
  rbind(gwas_out$Female_forsmani$GWAResult$phenos_group %>% 
          filter(!is.na(pValue)) %>% 
          mutate(Sex = "Female", Species = "forsmani")) %>% 
  rbind(gwas_out$Female_praehirsuta$GWAResult$phenos_group %>% 
          filter(!is.na(pValue)) %>% 
          mutate(Sex = "Female", Species = "praehirsuta")) %>% 
  rename(Chromosome = chr,
         Position = pos,
         Trait = trait) %>% 
  mutate(log_pval = -log10(pValue))
  


(GWAS_out %>% 
    filter(Sex == "Male") %>% 
    geom_manhattan(aes(y = log_pval, facetting1 = Species, facetting2 = Trait), alpha = 0.5)) +
  facet_grid2(Species ~ Trait, scales = "free")


