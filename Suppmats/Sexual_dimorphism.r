#' This script is written to look at the size distribution of the adults in the
#' experiment


########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "ggtext")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################################
################### Useful variables ###################
########################################################
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20),
        legend.text = element_markdown(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


########################################################
##################### Import data ######################
########################################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "parents") %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b"),
         animal = as.integer(str_remove_all(ID_DNA_RAD, "E"))) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()

########################################################
################## Plot distributions ##################
########################################################
size_distribution <- metadata %>% 
  mutate(Sex = ifelse(Sex == "F", "Female", "Male")) %>% 
  select(Sex, Species, Size) %>% 
  group_by(Sex, Species) %>% 
  mutate(group = row_number()) %>% 
  pivot_wider(names_from = Species, values_from = Size, names_prefix = "Size_") %>% 
  ungroup() %>% 
  select(-group) %>% 
  ggplot() +
  geom_histogram(data = . %>% 
                   mutate(Species = "forsmani"),
                 aes(x = Size_forsmani, y = ..count.., colour = Sex, fill = Species),
                 alpha = 0.7, position = "identity", lwd = 1) + 
  geom_histogram(data = . %>% 
                   mutate(Species = "praehirsuta"),
                 aes(x = Size_praehirsuta, y = -..count.., colour = Sex, fill = Species),
                 alpha = 0.7, position = "identity", lwd = 1) +
  scale_colour_manual(values = colours_sexes) +
  scale_fill_manual(values = colours_species,
                    labels = c("forsmani" = "*J. forsmani*",
                               "praehirsuta" = "*J. praehirsuta*")) +
  labs(x = "Size\n(mm)",
       y = "Count") +
  my_theme

size_distribution
ggsave(plot = size_distribution,
       filename = "../Figures/Suppfigs/Sexual_dimorphism_per_sp.png",
       width = 1000, height = 700, units = "px", scale = 4)


########################################################
################### Statistical test ###################
########################################################
aov(Size ~ Sex * Species, data = metadata) %>% 
  summary()
