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
  unique() %>% 
  ggplot() +
  geom_violin(aes(x = Species, y = Time_first_brood, fill = Species),
               alpha = 0.7, colour = "black",
              quantile.linetype = c(2, 1, 2)) +
  geom_point(aes(x = Species, y = Time_first_brood, fill = Species), pch = 21,
             size = 3, position = position_jitter(width = 0.25, height = 0),
             alpha = 0.7) +
  scale_fill_manual(values = colours_species,
                    labels = c("forsmani" = "J. forsmani",
                               "praehirsuta" = "J. praehirsuta")) +
  scale_x_discrete(labels = c("forsmani" = "J. forsmani",
                              "praehirsuta" = "J. praehirsuta")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y = "Time to produce first brood\n(days)") + 
  my_theme +
  theme(legend.text = element_text(face = "italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))

distribution_time_first_brood
ggsave(plot = distribution_time_first_brood,
       filename = "../Figures/Suppfigs/Distribution_time_first_brood.png",
       width = 1000, height = 700, scale = 4, units = "px")
