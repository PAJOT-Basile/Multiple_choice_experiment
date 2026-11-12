#' This script is written to look at the relationship in females between the 
#' size of the individual and the size of the brood it produced.


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

# Reproductive success
load("Associated_data/Rdata/Fitness_components_parents.rdata")
fitness_components <- fitness_components %>% 
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size),
            by = join_by("Parent_ID" == "ID_DNA_RAD")) %>% 
  group_by(Species, Sex) %>% 
  rename(W1 = nb_partners,
         W12 = nb_offspring,
         W123 = nb_offspring_surv60) %>% 
  select(-c(max_nb_offspring, nb_offspring_surv40, -contains("trunc"))) %>% 
  mutate(
    # Number of individuals
    p0 = 1/ n(),
    # Standardised fitnesses
    w1 = W1 / mean(W1, na.rm = TRUE),
    w12 = W12 / mean(W12, na.rm = TRUE),
    w123 = W123 / mean(W123, na.rm = TRUE),
    # Absolute fitnesses
    W2 = W12 / W1,
    W3 = W123 / W12,
    # REST
    p1 = p0 * w1,
    w2 = W2 / sum(p1 * W2, na.rm = TRUE),
    p2 = p1 * w2,
    w3 = W3 / sum(p2 * W3, na.rm = TRUE),
    p3 = p2 * w3,
    across(c(w1, w12, w123, p1, w2, p2, w3, p3), ~ ifelse(is.na(.), 0, .))) %>% 
  ungroup()

########################################################
################## Plot relationship ###################
########################################################
relationship <- fitness_components %>% 
  filter(Sex == "Female", W1 != 0) %>% 
  ggplot(aes(x = Size, y = W12, colour = Species)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(se = FALSE, method = "lm", lwd = 1.1) +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  labs(x = "Female size\n(mm)",
       y = "Number of offspring\n(W12)") + 
  my_theme +
  theme(legend.text = element_text(face = "italic"))

relationship
ggsave(plot = relationship,
       filename = "../Figures/Suppfigs/Female_size_brood_size.png",
       scale = 4, width = 800, height = 600, units = "px")
