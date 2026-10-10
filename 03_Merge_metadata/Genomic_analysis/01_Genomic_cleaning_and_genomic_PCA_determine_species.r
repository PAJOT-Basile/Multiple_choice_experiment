#########################################################
####################### Libraries #######################
#########################################################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet", "patchwork", "ggforce", "ggnewscale")
pacman::p_load(char = libraries, character.only = TRUE)

#########################################################
#################### Useful functions ###################
#########################################################
source("98_Transversal_scripts/Transversal_functions.r")

#########################################################
#################### Useful variables ###################
#########################################################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

colours_species <- c("praehirsuta" = "navy",
                    "forsmani" = "#3A9AB2",
                    "Intermediate" = "red")

set.seed(123456789)

#########################################################
###################### Import data ######################
#########################################################
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
         Phenotype = factor(Phenotype, levels = c("forsmani", "Intermediate", "praehirsuta")))

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

load("Associated_data/Rdata/Families_first_brood.rdata")
#########################################################
######################## Run PCA ########################
#########################################################
# The PCA will only be used to identify the species of the mothers (female parents)
# that were put in the aquarium. So, the PCA has to be run only on the parents
## Subset the genind to keep only the parents
parent_ids <- indNames(data) %>% 
  as.data.frame() %>% 
  rename(ID_DNA_RAD = ".") %>% 
  left_join(metadata, by = "ID_DNA_RAD") %>% 
  filter(Family_level == "parents") %>% 
  pull(ID_DNA_RAD)

parents_data <- data[indNames(data) %in% parent_ids]
pca <- scaleGen(parents_data, NA.method="mean",scale=F,center=T) %>%
  dudi.pca(scale=T, nf = 5,scannf = F)

# Get the variance explained by the different axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax4 <- ((pca$eig[4] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax5 <- ((pca$eig[5] / sum(pca$eig)) * 100) %>% round(digits = 2)


# Plot the results
pca_plot <- pca$li %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(metadata, by = "ID_DNA_RAD") %>% 
  mutate(Sex_pheno = case_when(Sex == "F" ~ "Females",
                               Phenotype == "forsmani" ~ "Males forsmani",
                               Phenotype == "praehirsuta" ~ "Males praehirsuta") %>%
           factor(levels = c("Males forsmani", "Males praehirsuta", "Females"))) %>%
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
  labs(x = paste0("Axis 1 (", var_ax1, "%)"),
       y = paste0("Axis 2 (", var_ax2, "%)")) +
  my_theme

pca_plot
ggsave("../Figures/Suppfigs/PCA_parents.png", pca_plot,
       scale = 4, width = 1200, height = 800, units = "px")

# The first axis of the PCA clearly separates the two species. To make sure that
# the identification of the species of female parents, we use the function find.clusters
# from the adegenet package allowing to identify genetic clusters on the PCA.
# We only need two clusters to differentiate the species, but we also try it out
# with four clusters to make sure the identified clusters are stable.

# Find clusters on the PCA
# For K = 2
clustering_pca2 <- find.clusters(parents_data, n.pca = 5, n.clust = 2)
# For K = 4
clustering_pca4 <- find.clusters(parents_data, n.pca = 5, n.clust = 4)

# Plot the results
# For K = 2
pca_2_groups <- clustering_pca2$grp %>% 
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(pca$li %>% 
              rownames_to_column("ID_DNA_RAD"),
            by = "ID_DNA_RAD") %>% 
  left_join(metadata,
            by = "ID_DNA_RAD") %>% 
  mutate(Group = factor(Group),
         Sex_pheno = case_when(Sex == "F" ~ "Females",
                               Phenotype == "forsmani" ~ "Males forsmani",
                               Phenotype == "praehirsuta" ~ "Males praehirsuta") %>% 
           factor(levels = c("Males forsmani", "Males praehirsuta", "Females"))) %>% 
  ggplot(aes(x = Axis1, y = -Axis2)) +
  geom_mark_ellipse(aes(colour = Group, label = Group), lwd = 1.2, label.fontsize = 15) +
  scale_colour_manual(values = c("2" = "navy", "1" = "#3A9AB2"),
                      guide = "none") +
  new_scale_colour() +
  geom_point(aes(colour = Sex_pheno, shape = Sex_pheno), size = 4, alpha = 0.7) +
  scale_colour_manual(name = "Sex and\nPhenotype",
                      values = c("Males forsmani" = "#3A9AB2",
                                 "Males praehirsuta" = "navy",
                                 "Females" = "grey50")) +
  scale_shape_manual(name = "Sex and\nPhenotype",
                     values = c("Males forsmani" = 17,
                                "Males praehirsuta" = 17,
                                "Females" = 15)) +
  labs(x = paste0("Axis 1 (", var_ax1, "%)"),
       y = paste0("Axis 2 (", var_ax2, "%)"),
       tag = "(A)") +
  xlim(-100, 100) +
  ylim(-95, 250) +
  my_theme +
  theme(legend.box.margin = margin(r = 2))

# For K = 4
pca_4_groups <- clustering_pca4$grp %>% 
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(pca$li %>% 
              rownames_to_column("ID_DNA_RAD"),
            by = "ID_DNA_RAD") %>% 
  left_join(metadata,
            by = "ID_DNA_RAD") %>% 
  mutate(Group = factor(Group),
         Sex_pheno = case_when(Sex == "F" ~ "Females",
                               Phenotype == "forsmani" ~ "Males forsmani",
                               Phenotype == "praehirsuta" ~ "Males praehirsuta") %>% 
           factor(levels = c("Males forsmani", "Males praehirsuta", "Females"))) %>%
  ggplot(aes(x = Axis1, y = -Axis2)) +
  geom_mark_ellipse(aes(colour = Group, label = Group), lwd = 1.2, label.fontsize = 15) +
  scale_colour_manual(values = c("2" = "navy", "1" = "#3A9AB2", "3" = "#7CAE00", "4" = "#C77CFF"),
                      guide = "none") +
  new_scale_colour() +
  geom_point(aes(colour = Sex_pheno, shape = Sex_pheno), size = 4, alpha = 0.7) +
  scale_colour_manual(name = "Sex and\nPhenotype",
                      values = c("Males forsmani" = "#3A9AB2",
                                 "Males praehirsuta" = "navy",
                                 "Females" = "grey50")) +
  scale_shape_manual(name = "Sex and\nPhenotype",
                     values = c("Males forsmani" = 17,
                                "Males praehirsuta" = 17,
                                "Females" = 15)) +
  labs(x = paste0("Axis 1 (", var_ax1, "%)"),
       tag = "(B)") +
  xlim(-100, 100) +
  ylim(-95, 250) +
  my_theme +
  theme(axis.title.y = element_blank(),
        legend.box.margin = margin(r = 2))


(pca_2_groups + pca_4_groups +
  plot_layout(guides = "collect", widths = c(1, 0.95)) &
  theme(legend.box.background = element_rect(colour = "white", size = 1))) %>% 
  ggsave(plot = ., filename = "../Figures/Suppfigs/PCA_parents_identification.png",
         width = 1150, height = 550, units = "px", scale = 4)


#########################################################
########### Get the species names from the PCA ##########
#########################################################
correspondence_species_group <- clustering_pca2$grp %>% 
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(metadata,
            by = "ID_DNA_RAD") %>%
  filter(Sex == "M") %>%
  group_by(Group) %>%
  summarize(Species = unique(Phenotype))

Species_names_parents <- clustering_pca2$grp %>% 
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(correspondence_species_group) %>% 
  select(ID_DNA_RAD, Species)

save(Species_names_parents, file = "Associated_data/Rdata/Correspondence_name_species.rdata")

#########################################################
##################### PCA on everyone ###################
#########################################################
pca <- scaleGen(data, NA.method = "mean", scale = FALSE, center = TRUE) %>%
  dudi.pca(scale = TRUE, nf = 5,scannf = FALSE)

# Get the variance explained by the different axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax4 <- ((pca$eig[4] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax5 <- ((pca$eig[5] / sum(pca$eig)) * 100) %>% round(digits = 2)

# Find clusters on the PCA
# For K = 2
clustering_pca2 <- find.clusters(data, n.pca = 5, n.clust = 2)



# Plot the results
pca_plot <- clustering_pca2$grp %>%
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(pca$li %>% 
              rownames_to_column("ID_DNA_RAD"),
            by = "ID_DNA_RAD") %>% 
  left_join(metadata,
            by = "ID_DNA_RAD") %>% 
  left_join(Families %>% 
              rbind(Species_names_parents),
            by = "ID_DNA_RAD") %>% 
  rename(Species_parents = Species) %>% 
  mutate(Group = factor(Group)) %>% 
  ggplot(aes(x = Axis1, y = Axis2)) +
  geom_mark_ellipse(aes(colour = Group, label = Group), lwd = 1.2, label.fontsize = 15) +
  scale_colour_manual(values = c("1" = "navy", "2" = "#3A9AB2"),
                      guide = "none") +
  new_scale_colour() +
  geom_point(aes(colour = Phenotype, shape = Species_parents), size = 4,
             alpha = 0.7) +
  scale_colour_manual(values = colours_species) +
  scale_shape_manual(name = "Parent's\nspecies",
                     values = c("forsmani" = 21, "praehirsuta" = 22)) +
  labs(x = paste0("Axis 1 (", var_ax1, "%)"),
       y = paste0("Axis 2 (", var_ax2, "%)")) +
  my_theme 

pca_plot
ggsave("../Figures/Suppfigs/PCA_everyone.png", pca_plot,
       scale = 3, width = 1200, height = 800, units = "px")


#########################################################
##################### PCA on offspring ##################
#########################################################
offspring_id <- metadata %>% 
  filter(Family_level == "offspring") %>% 
  pull(ID_DNA_RAD)
pca <- scaleGen(data[indNames(data) %in% offspring_id], NA.method = "mean", scale = FALSE, center = TRUE) %>%
  dudi.pca(scale = TRUE, nf = 5,scannf = FALSE)

# Get the variance explained by the different axis
var_ax1 <- ((pca$eig[1] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax2 <- ((pca$eig[2] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax3 <- ((pca$eig[3] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax4 <- ((pca$eig[4] / sum(pca$eig)) * 100) %>% round(digits = 2)
var_ax5 <- ((pca$eig[5] / sum(pca$eig)) * 100) %>% round(digits = 2)

# Find clusters on the PCA
# For K = 2
clustering_pca2 <- find.clusters(data[indNames(data) %in% offspring_id], n.pca = 5, n.clust = 2)



# Plot the results
pca_plot <- clustering_pca2$grp %>%
  as.data.frame() %>% 
  rename(Group = ".") %>% 
  rownames_to_column("ID_DNA_RAD") %>% 
  left_join(pca$li %>% 
              rownames_to_column("ID_DNA_RAD"),
            by = "ID_DNA_RAD") %>% 
  left_join(metadata,
            by = "ID_DNA_RAD") %>% 
  left_join(Families,
            by = "ID_DNA_RAD") %>% 
  rename(Species_parents = Species) %>% 
  mutate(Group = factor(Group)) %>% 
  ggplot(aes(x = -Axis1, y = Axis2)) +
  geom_mark_ellipse(aes(colour = Group, label = Group), lwd = 1.2, label.fontsize = 15) +
  scale_colour_manual(values = c("1" = "navy", "2" = "#3A9AB2"),
                      guide = "none") +
  new_scale_colour() +
  geom_point(aes(colour = Phenotype, shape = Species_parents), size = 4,
             alpha = 0.7) +
  scale_colour_manual(values = colours_species) +
  scale_shape_manual(name = "Parent's\nspecies",
                     values = c("forsmani" = 21, "praehirsuta" = 22)) +
  labs(x = paste0("Axis 1 (", var_ax1, "%)"),
       y = paste0("Axis 2 (", var_ax2, "%)")) +
  my_theme 

pca_plot
ggsave("../Figures/Suppfigs/PCA_offspring.png", pca_plot,
       scale = 3, width = 1200, height = 800, units = "px")

