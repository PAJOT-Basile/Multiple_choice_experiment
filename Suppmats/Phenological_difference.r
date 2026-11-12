#' This script is written to look at the distribution of the time needed to 
#' produce a brood depending on the species.


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

# Import the detail of the time needed to produce a brood
survival <- read.table("Associated_data/Raw_metadata/Survival_offspring.tsv",
                       sep = "\t", header = TRUE)

########################################################
################# Prepare data to plot #################
########################################################
date_data <- survival %>%
  left_join(metadata %>%
              select(Label, Species, Name_offsprings_brood1),
            by = join_by("Mother_ID" == "Label")) %>% 
  separate(Name_offsprings_brood1, c("min_plate_br1", "max_plate_br1"), sep = "-") %>% 
  mutate(across(contains("plate"), ~ case_when(is.na(.) ~ str_remove_all(min_plate_br1, "O") %>% as.numeric(),
                                               TRUE ~ str_remove_all(., "O") %>% as.numeric()))) %>% 
  filter(Plate_ID >= min_plate_br1, Plate_ID <= max_plate_br1) %>% 
  select(Mother_ID, Isolation_date, Species) %>% 
  mutate(Isolation_date = as.Date(Isolation_date, format = "%d/%m/%Y"),
         Time_first_brood = Isolation_date - as.Date("04/06/2024", format = "%d/%m/%Y"))

########################################################
###################### Plot data #######################
########################################################
distribution_time_first_brood <- date_data %>% 
  group_by(Species) %>% 
  mutate(n = row_number()) %>% 
  pivot_wider(names_from = Species, values_from = Time_first_brood, names_prefix = "Time_") %>% 
  ggplot() +
  geom_histogram(data = . %>% 
                   mutate(Species = "forsmani"),
                 aes(x = Time_forsmani, y = ..count.., fill = Species),
                 alpha = 0.7, position = "identity", colour = "black") + 
  geom_histogram(data = . %>% 
                   mutate(Species = "praehirsuta"),
                 aes(x = Time_praehirsuta, y = -..count.., fill = Species),
                 alpha = 0.7, position = "identity", colour = "black") + 
  geom_hline(yintercept = 0, colour = "black", lwd = 0.8) + 
  scale_fill_manual(values = colours_species,
                    labels = c("forsmani" = "J. forsmani",
                               "praehirsuta" = "J. praehirsuta")) +
  labs(x = "Time to produce first brood\n(days)",
       y = "Count") + 
  my_theme +
  theme(legend.text = element_text(face = "italic"))

distribution_time_first_brood
ggsave(plot = distribution_time_first_brood,
       filename = "../Figures/Suppfigs/Distribution_time_first_brood.png",
       width = 1000, height = 700, scale = 4, units = "px")
