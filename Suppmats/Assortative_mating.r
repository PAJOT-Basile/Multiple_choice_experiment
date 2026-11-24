#' This script is written to check for assortative mating (check if larger males
#' reproduce with larger females or not)
#' 
#' Basile Pajot
#' Last updated: 17/11/2026


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

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20))


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

# Get the names of the parents in each couple that produced a brood
load("Associated_data/Rdata/Families_first_brood.rdata")
Families_first_brood <- Families_first_brood %>%
  select(Dad, Mom) %>% 
  unique() %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size, Sum_curved_setaeP1P5_default),
            by = join_by("Dad" == "ID_DNA_RAD")) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size, Species),
            by = join_by("Mom" == "ID_DNA_RAD"),
            suffix = c("_dad", "_mom"))


########################################################
###################### Plot data #######################
########################################################
mate_choice_size <- Families_first_brood %>% 
  ggplot(aes(x = Size_dad, y = Size_mom, colour = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  labs(x = "Size fathers\n(mm)",
       y = "Size mothers\n(mm)") +
  facet_wrap(vars(Species), scales = "free", nrow = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


mate_choice_size
ggsave(mate_choice_size,
       filename = "../Figures/Suppfigs/Relationship_size_female_size_male.png",
       units = "px", width = 1200, height = 800, scale = 4)

mate_choice_setae <- Families_first_brood %>% 
  ggplot(aes(x = Sum_curved_setaeP1P5_default, y = Size_mom, colour = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  labs(x = "Number of sexual setae",
       y = "Size mothers\n(mm)") +
  facet_wrap(vars(Species), scales = "free", nrow = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


mate_choice_setae
ggsave(mate_choice_setae,
       filename = "../Figures/Suppfigs/Relationship_size_female_setae_male.png",
       units = "px", width = 1200, height = 800, scale = 4)

####################################################################
###################### Plot size differences #######################
####################################################################


size_differences_realised_potential_couples <- expand_grid(metadata %>% 
              filter(Species == "forsmani", Sex == "F") %>% 
              pull(Size),
            metadata %>% 
              filter(Species == "forsmani", Sex == "M") %>% 
              pull(Size)) %>% 
  rename(Size_mom = "... %>% pull(Size)...1",
         Size_dad = "... %>% pull(Size)...2") %>% 
  mutate(Species = "forsmani") %>% 
  rbind(expand_grid(metadata %>% 
                      filter(Species == "praehirsuta", Sex == "F") %>% 
                      pull(Size),
                    metadata %>% 
                      filter(Species == "praehirsuta", Sex == "M") %>% 
                      pull(Size)) %>% 
          rename(Size_mom = "... %>% pull(Size)...1",
                 Size_dad = "... %>% pull(Size)...2") %>% 
          mutate(Species = "praehirsuta")) %>% 
  mutate(Group = "All possible") %>% 
  rbind(Families_first_brood %>% 
          select(Size_mom, Size_dad, Species) %>% 
          mutate(Group = "Observed")) %>% 
  mutate(Delta_size = Size_mom - Size_dad,
         Sp_group = case_when(Group == "All possible" ~ "All possible",
                              TRUE ~ paste(Group, Species))) %>% 
  ggplot(aes(x = Group, y = Delta_size, fill = Sp_group)) +
  geom_violin(quantile.linetype = c(2, 1, 2), alpha = 0.7, lwd = 1) +
  geom_point(data = . %>% 
               filter(Group != "All possible"), pch = 21,
             size = 3, position = position_jitter(width = 0.25, height = 0)) +
  scale_fill_manual(name = "Species",
                    values = c("All possible" = "grey40",
                               "Observed praehirsuta" = "navy",
                               "Observed forsmani" = "#3A9AB2"),
                    labels = c("All possible" = "All possible",
                               "Observed praehirsuta" = expression("Observed " * italic("J. praehirsuta")),
                               "Observed forsmani" = expression("Observed " * italic("J. forsmani"))),
                    drop = TRUE) +
  labs(y = "Size female - Size male\n(mm)") +
  facet_wrap(vars(Species), nrow = 1,
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  my_theme +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = "italic"))

size_differences_realised_potential_couples
ggsave("../Figures/Suppfigs/Size_differences_realised_potential_couples.png",
       plot = size_differences_realised_potential_couples,
       units = "px", scale = 4, width = 1000, height = 750)

