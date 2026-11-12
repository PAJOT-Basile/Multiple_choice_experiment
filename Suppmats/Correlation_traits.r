#' This script is written to look at the correlation of size with the other traits
#' studied here (the total number of spines and the total number of sexual setae)


########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse")
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
size_trait_relationship <- metadata %>% 
  filter(Sex == "M") %>% 
  select(ID_DNA_RAD, Species, Sex, Size, Sum_curved_setaeP1P5_default, Sum_spinesP6P7) %>% 
  pivot_longer(contains("Sum"), names_to = "Trait2", values_to = "Count") %>% 
  ggplot(aes(x = Size, y = Count, colour = Species)) + 
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.1) +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  facet_wrap(vars(Trait2), scales = "free_y",
             labeller = as_labeller(c(
               "Sum_curved_setaeP1P5_default" = "Number of sexual setae",
               "Sum_spinesP6P7" = "Number of spines"
             )),
            strip.position = "left",
            ncol = 1) +
  labs(x = "Size\n(mm)") +
  my_theme +
  theme(legend.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, colour = NA))

size_trait_relationship
ggsave(plot = size_trait_relationship,
       filename = "../Figures/Suppfigs/Size_trait_covariance.png",
       width = 800, height = 1000, units = "px", scale = 4)
