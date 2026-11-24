#' This script is written to plot parent-offspring regressions to look at the 
#' raw data.
#' 
#' Basile Pajot
#' Last updated: 23/11/2026


########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################
########### Useful variables ###########
########################################
set.seed(123456789)
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20))

########################################
############# Import data ##############
########################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b"),
         animal = as.integer(str_remove_all(ID_DNA_RAD, "E"))) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(  # For the offspring, the number of setae was counted more precisely than for the parents, 
    # especially on the carps of the different pereiopods. For the values to be 
    # comparable between parents and offspring, we remove the carps from calculations.
    Sum_curved_setaeP1P5_trunc = Sum_curved_setaeP1P5_default - default_P1.carp.curv.setae - default_P2.carp.curv.setae - default_P3.carp.curv.setae - default_P4.carp.curv.setae - default_P5.carp.curv.setae,
    across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup() %>% 
  filter(!is.na(Sex))

# Make the pedigree
load("Associated_data/Rdata/Families_first_brood.rdata")
pedigree_with_phenotypes <- Families_first_brood %>%
  as_tibble() %>% 
  select(-c(nb_tot_juvs, Brood_name)) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Sex, Size, Sum_curved_setaeP1P5_trunc, Sum_spinesP6P7),
            by = join_by("Offspring" == "ID_DNA_RAD")) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size, Sum_curved_setaeP1P5_trunc, Sum_spinesP6P7),
            by = join_by("Dad" == "ID_DNA_RAD"),
            suffix = c("_offspring", "_dad")) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size, Species),
            by = join_by("Mom" == "ID_DNA_RAD")) %>% 
  rename(Size_mom = Size)


size_table <- pedigree_with_phenotypes %>% 
  group_by(Mom, Dad) %>% 
  summarize(Size_mom = unique(Size_mom),
            Size_dad = unique(Size_dad),
            Mean_size_offspring = mean(Size_offspring, na.rm = TRUE),
            Species = unique(Species)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(Mean_size_parents = mean(c(Size_mom, Size_dad)))

setae_spines_table <- pedigree_with_phenotypes %>% 
  group_by(Mom, Dad) %>% 
  summarize(Setae_dad = unique(Sum_curved_setaeP1P5_trunc_dad),
            Spines_dad = unique(Sum_spinesP6P7_dad),
            Mean_setae_offspring = mean(Sum_curved_setaeP1P5_trunc_offspring, na.rm = TRUE),
            Mean_spines_offspring = mean(Sum_spinesP6P7_offspring, na.rm = TRUE),
            Species = unique(Species)) %>% 
  ungroup()


##############################################
############## Plot data setae ###############
##############################################
min_nb_setae_prae <- setae_spines_table %>% 
  filter(Species == "praehirsuta") %>% 
  select(Mean_setae_offspring, Setae_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Setae") %>% 
  pull(Setae) %>% 
  min(na.rm = TRUE)
max_nb_setae_prae <- setae_spines_table %>% 
  filter(Species == "praehirsuta") %>% 
  select(Mean_setae_offspring, Setae_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Setae") %>% 
  pull(Setae) %>% 
  max(na.rm = TRUE)

min_nb_setae_forsmani <- setae_spines_table  %>%
  filter(Species == "forsmani") %>%  
  select(Mean_setae_offspring, Setae_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Setae") %>% 
  pull(Setae) %>% 
  min(na.rm = TRUE)
max_nb_setae_forsmani <- setae_spines_table  %>%
  filter(Species == "forsmani") %>%  
  select(Mean_setae_offspring, Setae_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Setae") %>% 
  pull(Setae) %>% 
  max(na.rm = TRUE)


parent_offspring_regression_setae <- setae_spines_table %>%
  ggplot(aes(x = Setae_dad, y = Mean_setae_offspring, colour = Species)) +
  geom_line(aes(group = Mom), colour = "black", alpha = 0.4,
            lwd = 0.9) +
  geom_line(aes(group = Dad), colour = "red", alpha = 0.4,
            lwd = 0.9) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_smooth(data = data.frame(x = c(seq(min_nb_setae_prae, max_nb_setae_prae, 0.01),
                                      seq(min_nb_setae_forsmani, max_nb_setae_forsmani, 0.01)),
                                y = c(seq(min_nb_setae_prae, max_nb_setae_prae, 0.01),
                                      seq(min_nb_setae_forsmani, max_nb_setae_forsmani, 0.01)),
                                Species = c(rep("praehirsuta", length(seq(min_nb_setae_prae, max_nb_setae_prae, 0.01))),
                                            rep("forsmani", length(seq(min_nb_setae_forsmani, max_nb_setae_forsmani, 0.01))))),
              aes(x = x, y = y), colour = "black", se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  # coord_cartesian(xlim = c(min_nb_setae, max_nb_setae),
  #                 ylim = c(min_nb_setae, max_nb_setae)) +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Number of sexual setae parent",
       y = "Mean number of sexual setae for offspring") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


parent_offspring_regression_setae
ggsave("../Figures/Suppfigs/Parent_offspring_regression_setae.png",
       parent_offspring_regression_setae, units = "px", scale = 4,
       height = 500, width = 1000)



###############################################
############## Plot data spines ###############
###############################################
min_nb_spines <- setae_spines_table %>% 
  select(Mean_spines_offspring, Spines_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Spines") %>% 
  pull(Spines) %>% 
  min(na.rm = TRUE)
max_nb_spines <- setae_spines_table %>% 
  select(Mean_spines_offspring, Spines_dad) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Spines") %>% 
  pull(Spines) %>% 
  max(na.rm = TRUE)


parent_offspring_regression_spines <- setae_spines_table %>%
  ggplot(aes(x = Spines_dad, y = Mean_spines_offspring, colour = Species)) +
  geom_line(aes(group = Mom), colour = "black", alpha = 0.4,
            lwd = 0.9) +
  geom_line(aes(group = Dad), colour = "red", alpha = 0.4,
            lwd = 0.9) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_smooth(data = data.frame(x = seq(min_nb_spines, max_nb_spines, 0.01),
                                y = seq(min_nb_spines, max_nb_spines, 0.01)),
              aes(x = x, y = y), colour = "black", se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  coord_cartesian(xlim = c(min_nb_spines, max_nb_spines),
                  ylim = c(min_nb_spines, max_nb_spines)) +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Number of spines parent",
       y = "Mean number of spines for offspring") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


parent_offspring_regression_spines
ggsave("../Figures/Suppfigs/Parent_offspring_regression_spines.png",
       parent_offspring_regression_spines, units = "px", scale = 4,
       height = 500, width = 1000)

#############################################
############## Plot data size ###############
#############################################
min_mean_size <- size_table %>% 
  select(Mean_size_parents, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  min(na.rm = TRUE)

max_mean_size <- size_table %>% 
  select(Mean_size_parents, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  max(na.rm = TRUE)

parent_offspring_regression_size <- size_table %>%
  ungroup() %>% 
  ggplot(aes(x = Mean_size_parents, y = Mean_size_offspring, colour = Species)) +
  geom_line(aes(group = Mom), colour = "black", alpha = 0.4,
            lwd = 0.9) +
  geom_line(aes(group = Dad), colour = "red", alpha = 0.4,
            lwd = 0.9) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_smooth(data = data.frame(x = seq(min_mean_size, max_mean_size, 0.01),
                                y = seq(min_mean_size, max_mean_size, 0.01)),
              aes(x = x, y = y), colour = "black", se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  coord_cartesian(xlim = c(min_mean_size, max_mean_size),
                  ylim = c(min_mean_size, max_mean_size)) +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Mean parent size",
       y = "Mean offspring size") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


parent_offspring_regression_size
ggsave("../Figures/Suppfigs/Parent_offspring_regression_size.png",
       parent_offspring_regression_size, units = "px", scale = 4,
       height = 500, width = 1000)

##################################################################
############## Plot data size with separated sexes ###############
##################################################################
size_table_sexes <- pedigree_with_phenotypes %>% 
  group_by(Mom, Dad, Sex) %>% 
  summarize(Size_mom = unique(Size_mom),
            Size_dad = unique(Size_dad),
            Mean_size_offspring = mean(Size_offspring, na.rm = TRUE),
            Species = unique(Species)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(Mean_size_parent = case_when(Sex == "F" ~ Size_mom,
                                      TRUE ~ Size_dad))


min_mean_size_males <- size_table_sexes %>% 
  filter(Sex == "M") %>% 
  select(Mean_size_parent, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  min(na.rm = TRUE)

max_mean_size_males <- size_table_sexes %>% 
  filter(Sex == "M") %>% 
  select(Mean_size_parent, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  max(na.rm = TRUE)

parent_offspring_regression_size_males <- size_table_sexes %>%
  ungroup() %>% 
  filter(Sex == "M") %>% 
  ggplot(aes(x = Mean_size_parent, y = Mean_size_offspring, colour = Species)) +
  geom_line(aes(group = Mom), colour = "black", alpha = 0.4,
            lwd = 0.9) +
  geom_line(aes(group = Dad), colour = "red", alpha = 0.4,
            lwd = 0.9) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_smooth(data = data.frame(x = seq(min_mean_size_males, max_mean_size_males, 0.01),
                                y = seq(min_mean_size_males, max_mean_size_males, 0.01)),
              aes(x = x, y = y), colour = "black", se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  coord_cartesian(xlim = c(min_mean_size_males, max_mean_size_males),
                  ylim = c(min_mean_size_males, max_mean_size_males)) +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Mean father size",
       y = "Mean male offspring size") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


parent_offspring_regression_size_males
ggsave("../Figures/Suppfigs/Parent_offspring_regression_size_only_males.png",
       parent_offspring_regression_size_males, units = "px", scale = 4,
       height = 500, width = 1000)


min_mean_size_females <- size_table_sexes %>% 
  filter(Sex == "F") %>% 
  select(Mean_size_parent, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  min(na.rm = TRUE)

max_mean_size_females <- size_table_sexes %>% 
  filter(Sex == "F") %>% 
  select(Mean_size_parent, Mean_size_offspring) %>% 
  drop_na() %>% 
  pivot_longer(everything(), names_to = "Family_level", values_to = "Size") %>% 
  pull(Size) %>% 
  max(na.rm = TRUE)

parent_offspring_regression_size_females <- size_table_sexes %>%
  ungroup() %>% 
  filter(Sex == "F") %>% 
  ggplot(aes(x = Mean_size_parent, y = Mean_size_offspring, colour = Species)) +
  geom_line(aes(group = Mom), colour = "black", alpha = 0.4,
            lwd = 0.9) +
  geom_line(aes(group = Dad), colour = "red", alpha = 0.4,
            lwd = 0.9) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_smooth(data = data.frame(x = seq(min_mean_size_females, max_mean_size_females, 0.01),
                                y = seq(min_mean_size_females, max_mean_size_females, 0.01)),
              aes(x = x, y = y), colour = "black", se = FALSE, method = "lm") +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta")) +
  coord_cartesian(xlim = c(min_mean_size_females, max_mean_size_females),
                  ylim = c(min_mean_size_females, max_mean_size_females)) +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  labs(x = "Mean mother size",
       y = "Mean female offspring size") +
  my_theme +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))


parent_offspring_regression_size_females
ggsave("../Figures/Suppfigs/Parent_offspring_regression_size_only_females.png",
       parent_offspring_regression_size_females, units = "px", scale = 4,
       height = 500, width = 1000)

