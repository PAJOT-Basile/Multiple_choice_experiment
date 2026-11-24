#' This script is written to look at the number of spines of different sizes and 
#' the correlation with the size of the males.
#' 
#' Basile Pajot
#' Last updated: 17/11/2026


########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "ggh4x")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################################
################### Useful variables ###################
########################################################
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20),
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
###################### Plot data #######################
########################################################
size_spine_size_indiv <- metadata %>% 
  filter(Sex == "M") %>% 
  mutate(Large_spines = P7.larg.ep + P6.larg.ep,
         Medium_spines = P7.med.ep + P6.med.ep,
         Small_spines = P7.small.ep + P6.small.ep) %>% 
  select(Species, Size, Large_spines, Medium_spines, Small_spines) %>% 
  pivot_longer(contains("spines"), names_to = "Spine_size", values_to = "Count") %>% 
  ggplot(aes(x = Size, y = Count, colour = Species)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  labs(x = "Size\n(mm)") +
  facet_grid2(Spine_size ~ Species,
             labeller = as_labeller(c(
               "Large_spines" = "Number of\nlarge spines",
               "Medium_spines" = "Number of\nmedium spines",
               "Small_spines" = "Number of\nsmall spines",
               "forsmani" = "J. forsmani",
               "praehirsuta" = "J. praehirsuta")),
             switch = "y") +
  my_theme +
  theme(legend.text = element_text(face = "italic"),
        strip.text.x = element_text(face = "italic"),
        axis.title.y = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 20),
        strip.placement = "outside")

size_spine_size_indiv
ggsave("../Figures/Suppfigs/Relationship_size_indiv_size_spines.png",
       size_spine_size_indiv, scale = 4, units = "px",
       height = 800, width = 1000)
