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
libraries <- c("tidyverse", "adegenet", "vcfR", "statgenGWAS", "ggnewscale", "patchwork")
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

vcf_file <- "Associated_data/VCF_File.vcf"
# vcftools <- "/shared/software/miniconda/envs/vcftools-0.1.16/bin/vcftools"

group_names <- c("Male_forsmani", "Male_praehirsuta", "Female_forsmani", "Female_praehirsuta")

##############################
###### Useful functions ######
##############################
source("98_Transversal_scripts/Transversal_functions.r")

##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "offspring") %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b"),
         animal = as.integer(str_remove_all(ID_DNA_RAD, "E"))) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
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
# metadata %>% 
#   filter(Sex == "M", Species == "forsmani") %>% 
#   ggplot(aes(x = Size, y = Sum_spinesP6P7)) +
#   geom_smooth(method = "lm", se = FALSE, formula = y ~ x + I(x^2)) +
#   geom_point()
# # Males forsmani spines
# lm_forsmani_spines <- metadata %>%
#   filter(Sex == "M", Species == "forsmani") %>% 
#   lm(Sum_spinesP6P7 ~ Size + I(Size^2), data = .)
# 
# # Males forsmani setae
# lm_forsmani_setae <- metadata %>%
#   filter(Sex == "M", Species == "forsmani") %>% 
#   lm(Sum_curved_setaeP1P5_default  ~ Size, data = .)
# 
# # Males praehirsuta setae
# lm_praehirsuta_setae <- metadata %>%
#   filter(Sex == "M", Species == "praehirsuta") %>% 
#   lm(Sum_curved_setaeP1P5_default ~ Size, data = .)
#   
# # Males praehirsuta spines
# lm_praehirsuta_spines <- metadata %>% 
#   filter(Sex == "M", Species == "praehirsuta") %>% 
#   lm(Sum_spinesP6P7 ~ Size, data = .)
# 
# lms <- list("forsmani_spines" = lm_forsmani_spines,
#             "forsmani_setae" = lm_forsmani_setae,
#             "praehirsuta_spines" = lm_praehirsuta_spines,
#             "praehirsuta_setae" = lm_praehirsuta_setae)


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
                  select(ID_DNA_RAD, Sex, Size, Sum_curved_setaeP1P5_default, Sum_spinesP6P7),
                by = join_by("INDV1" == "ID_DNA_RAD")) %>% 
      left_join(metadata %>% 
                  select(ID_DNA_RAD, Sex, Size, Sum_curved_setaeP1P5_default, Sum_spinesP6P7),
                by = join_by("INDV2" == "ID_DNA_RAD")) %>% 
      filter(Sex.x == str_split(sex, "", simplify = TRUE)[, 1],
             Sex.y == str_split(sex, "", simplify = TRUE)[, 1],
             !is.na(Size.x), !is.na(Size.y),
             !is.na(Sum_curved_setaeP1P5_default.x), !is.na(Sum_curved_setaeP1P5_default.y),
             !is.na(Sum_spinesP6P7.x), !is.na(Sum_spinesP6P7.y)) %>% 
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
      select(ID_DNA_RAD, Size, Sum_curved_setaeP1P5_default, Sum_spinesP6P7) %>% 
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
gwas_out <- lapply(group_names, function(group_name, geno, map_chromosome, phenos, scaled_relatedness){
  print(group_name)
  if (grepl("Female", group_name)){
    phenos_group <- phenos[[group_name]] %>% 
      select(-Sum_curved_setaeP1P5_default)
  }else{
    phenos_group <- phenos[[group_name]]
  }
  gdata_object <- createGData(geno = geno,
              map = map_chromosome,
              pheno = phenos_group,
              kin = scaled_relatedness[[group_name]])

  runSingleTraitGwas(gdata_object, thrType = "bonf")
  
}, geno, map_chromosome, phenos, scaled_relatedness)
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
  

#####################################################
############## Prepare Manhattan plot ###############
#####################################################
# TODO: tout préparer pour tester de faire un joli manhattan plot
nb_chromosomes <- GWAS_out$Chromosome %>%
  unique() %>%
  length()


colour_chromosome <- data.frame(
  "Chromosome" = rep(GWAS_out$Chromosome %>% unique(), 2),
  "Colour" = c(rep(c("#3A9AB2", "#9EE1EB"), nb_chromosomes)[1:nb_chromosomes],
               rep(c("navy", "#A0A0D7"), nb_chromosomes)[1:nb_chromosomes]),
  "Species" = c(rep("forsmani", nb_chromosomes),
                rep("praehirsuta", nb_chromosomes))
)

default_colours <- data.frame(
  "Chromosome" = GWAS_out$Chromosome %>% unique(),
  "Colour" = rep(c("grey41", "grey4"), nb_chromosomes)[1:nb_chromosomes]
)
  

GWAS_out_prep <- GWAS_out %>% 
  mutate(Locus_name = str_remove_all(snp, ":[+-].[01]$")) %>% 
  left_join(summary_stats %>% 
              select(Locus_name, BP_cumul),
            by = "Locus_name",
            relationship = "many-to-one")


chromosome_centers <- GWAS_out_prep %>% 
  filter(Species == "forsmani") %>% 
  group_by(Chromosome) %>% 
  # The values of the centers of the chromosomes are approximated using the
  # mean function
  summarise(center = mean(BP_cumul, na.rm = TRUE))



threshold_pvalue <- 5

##################################################
############## Plot primary figure ###############
##################################################
GWAS_setae_plot_forsmani <- GWAS_out_prep %>%
  filter(Sex == "Male",
         Trait == "Sum_curved_setaeP1P5_default",
         Species == "forsmani") %>% 
  ggplot() +
  geom_point(data = . %>% 
               filter(log_pval < threshold_pvalue), size = 3, alpha = 0.2,
             aes(x = BP_cumul, y = log_pval, colour = Chromosome)) +
  scale_colour_manual(values = default_colours %>% 
                        pull(Colour),
                      guide = "none") +
  new_scale_colour() +
  geom_point(data = . %>% 
               filter(log_pval > threshold_pvalue), size = 3, alpha = 0.9,
             aes(x = BP_cumul, y = log_pval, colour = Chromosome)) +
  scale_colour_manual(values = colour_chromosome %>% 
                        filter(Species == "forsmani") %>% 
                        pull(Colour),
                      guide = "none") +
  geom_hline(yintercept = threshold_pvalue, colour = "grey40", lty = "dotted",
             lwd = 1) +
  scale_x_continuous(labels = (chromosome_centers$Chromosome %>%
                                 str_split_fixed(., "G", 2))[, 2],
                     breaks = chromosome_centers$center) +
  labs(tag = "(A)") +
  facet_wrap(vars(Species),
             labeller = as_labeller(c("forsmani" = "J. forsmani"))) +
  my_theme +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        plot.tag = element_text(size = 20))

GWAS_setae_plot_praehirsuta <- GWAS_out_prep %>%
  filter(Sex == "Male",
         Trait == "Sum_curved_setaeP1P5_default",
         Species == "praehirsuta") %>% 
  ggplot() +
  geom_point(data = . %>% 
               filter(log_pval < threshold_pvalue), size = 3, alpha = 0.2,
             aes(x = BP_cumul, y = log_pval, colour = Chromosome)) +
  scale_colour_manual(values = default_colours %>% 
                        pull(Colour),
                      guide = "none") +
  new_scale_colour() +
  geom_point(data = . %>% 
               filter(log_pval > threshold_pvalue), size = 3, alpha = 0.9,
             aes(x = BP_cumul, y = log_pval, colour = Chromosome)) +
  scale_colour_manual(values = colour_chromosome %>% 
                        filter(Species == "praehirsuta") %>% 
                        pull(Colour),
                      guide = "none") +
  geom_hline(yintercept = threshold_pvalue, colour = "grey40", lty = "dotted",
             lwd = 1) +
  scale_x_continuous(labels = (chromosome_centers$Chromosome %>%
                                 str_split_fixed(., "G", 2))[, 2],
                     breaks = chromosome_centers$center) +
  labs(tag = "(C)") +
  facet_wrap(vars(Species),
             labeller = as_labeller(c("praehirsuta" = "J. praehirsuta"))) +
  my_theme +
  labs(x = "Chromosome") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        plot.tag = element_text(size = 20))



distribution_setae_forsmani <- metadata %>%
  filter(Species == "forsmani",
         Sex == "M") %>% 
  ggplot(aes(x = Sum_curved_setaeP1P5_default)) +
  geom_histogram(fill = "#3A9AB2", bins = 50) +
  facet_wrap(vars(Species),
             labeller = as_labeller(c("forsmani" = "J. forsmani"))) +
  labs(tag = "(B)") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        axis.title = element_blank(),
        plot.tag = element_text(size = 20))

distribution_setae_praehirsuta <- metadata %>%
  filter(Species == "praehirsuta",
         Sex == "M") %>% 
  ggplot(aes(x = Sum_curved_setaeP1P5_default)) +
  geom_histogram(fill = "navy", bins = 50) +
  facet_wrap(vars(Species),
             labeller = as_labeller(c("praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Sum of curved setae",
       tag = "(D)") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        plot.tag = element_text(size = 20))

title_y_axis <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold('-lo'*g[10] *' p-value for number of sexual setae')"),
            parse = TRUE, size = 6, angle = 90) +
  theme_void()


GWAS_setae_plot <- title_y_axis + GWAS_setae_plot_forsmani + distribution_setae_forsmani +
  GWAS_setae_plot_praehirsuta + distribution_setae_praehirsuta +
  plot_layout(design = "ABC
              ADE",
              widths = c(0.1, 1, 0.3))

GWAS_setae_plot

ggsave("../Figures/tmp/GWAS_setae.png", GWAS_setae_plot,
       scale = 4, width = 1200, height = 800, units = "px")


#########################################################
############## Plot supplementary figures ###############
#########################################################
GWAS_size_plot_m <- (GWAS_out %>% 
                     filter(Sex == "Male",
                            Trait == "Size") %>% 
                     geom_manhattan(aes(y = log_pval, facetting1 = Species), alpha = 0.5)) +
  facet_wrap(vars(Species), scales = "free", ncol = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  theme(strip.text = element_text(face = "italic"))

GWAS_size_plot_m
ggsave("../Figures/tmp/GWAS_size_males.png", GWAS_size_plot_m,
       scale = 4, width = 1200, height = 800, units = "px")


GWAS_size_plot_f <- (GWAS_out %>% 
                       filter(Sex == "Female",
                              Trait == "Size") %>% 
                       geom_manhattan(aes(y = log_pval, facetting1 = Species), alpha = 0.5)) +
  facet_wrap(vars(Species), scales = "free", ncol = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  theme(strip.text = element_text(face = "italic"))

GWAS_size_plot_f
ggsave("../Figures/tmp/GWAS_size_females.png", GWAS_size_plot_f,
       scale = 4, width = 1200, height = 800, units = "px")



GWAS_spines_plot <- (GWAS_out %>% 
                      filter(Sex == "Male",
                             Trait == "Sum_spinesP6P7") %>% 
                      geom_manhattan(aes(y = log_pval, facetting1 = Species), alpha = 0.5)) +
  facet_wrap(vars(Species), scales = "free", ncol = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  theme(strip.text = element_text(face = "italic"))

GWAS_spines_plot
ggsave("../Figures/tmp/GWAS_spines.png", GWAS_spines_plot,
       scale = 4, width = 1200, height = 800, units = "px")


