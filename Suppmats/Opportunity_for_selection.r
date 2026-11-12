# Figure to represent the opportunity for selection across the life cycle of the individuals


########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################
########### Useful variables ###########
########################################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())


colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("forsmani" = "#3A9AB2", "praehirsuta" = "navy")

########################################
########### Useful functions ###########
########################################
source("98_Transversal_scripts/Transversal_functions.r")

########################################
############# Import data ##############
########################################
# Import metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the dulicates
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  rowwise() %>% 
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all("a|b"),
         Family_level = factor(Family_level, levels = c("parents", "offspring")),
         Species = factor(Species, levels = c("forsmani", "hybrid", "praehirsuta"))) %>% 
  unique()


# Reproductive success
load("Associated_data/Rdata/Fitness_components_parents.rdata")
fitness_components <- fitness_components %>% 
  left_join(metadata %>%
              select(ID_DNA_RAD, Size),
            by = join_by("Parent_ID" == "ID_DNA_RAD")) %>% 
  # Here, we remove the female that had one offspring that was not genotyped
  filter(Parent_ID != "E4961") %>% 
  drop_na() %>% 
  group_by(Species, Sex) %>% 
  rename(W1 = nb_partners,
         W12 = nb_offspring,
         W123 = nb_offspring_surv60) %>% 
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
    across(c(w1, w12, w123, p1, w2, p2, w3, p3), ~ ifelse(is.na(.), 0, .)))

load("Associated_data/Rdata/distribution_variances_50000_boots_three_episodes.rda")


CI_stats <- distrib_opportunity_selection %>% 
  pivot_longer(cols = !c(Sex, Species), names_to = "Stat_name", values_to = "Stat") %>% 
  group_by(Species, Sex, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat, probs = 0.025),
            CI_97.5 = quantile(Stat, probs = 0.975))




########################################
## Compute opportunity for selection ###
########################################
I_vals <- fitness_components %>% 
  group_by(Species, Sex) %>% 
  summarize(I1 = sum(p0 * (w1 - 1)^2, na.rm = TRUE),
            I2 = sum(p1 * (w2 - 1)^2, na.rm = TRUE),
            I12 = sum(p0 * (w12 - 1)^2, na.rm = TRUE),
            
            COI1.2 = sum(p0 * (w1 - 1) * (w2 - 1), na.rm = TRUE),
            COI1.2_1 = sum(p1 * (w1 - 1) * (w2 - 1), na.rm = TRUE),
            COI12.2_1 = sum(p1 * (w12 - 1) * (w2 - 1), na.rm = TRUE),
            COI12.2 = sum(p0 * (w12 - 1) * (w2 - 1), na.rm = TRUE),
            
            I3 = sum(p2 * (w3 - 1)^2, na.rm = TRUE),
            I123 = sum(p0 * (w123 - 1)^2, na.rm = TRUE),
            
            COI1.3 = sum(p0 * (w1 - 1) * (w3 - 1), na.rm = TRUE),
            COI2.3 = sum(p1 * (w2 - 1) * (w3 - 1), na.rm = TRUE),
            COI12.3 = sum(p0 * (w12 - 1) * (w3 - 1), na.rm = TRUE),
            COI12.3_2 = sum(p2 * (w12 - 1) * (w3 - 1), na.rm = TRUE),
            COI123.3_2 = sum(p2 * (w123 - 1) * (w3 - 1), na.rm = TRUE),
            COI123.3 = sum(p0 * (w123 - 1) * (w3 - 1), na.rm = TRUE),
  ) %>% 
  mutate(C1 = COI1.2 + COI1.2_1 + (COI12.2_1 - COI12.2),
         C2 = COI12.3 + COI12.3_2 + (COI123.3_2 - COI123.3),
         across(c(contains("I", ignore.case = FALSE), contains("C")),
                ~ (. / I123) * 100,
                .names = "Percentage_{col}"),
         across(c(contains("I", ignore.case = FALSE), contains("C")), ~ round(., digits = 2))) %>% 
  select(-c(starts_with("Percentage"), starts_with("COI"))) %>%
  mutate(I_tot = I123) %>%
  pivot_longer(!c(Species, Sex, I_tot),
               names_to = "Stat_name",
               values_to = "Opportunity") %>%
  left_join(CI_stats,
            by = c("Species", "Sex", "Stat_name")) %>% 
  mutate(Stat_name = case_when(Stat_name == "I1" ~ "1\nNumber of mates",
                               Stat_name == "I2" ~ "2\nNumber of offspring per mate",
                               Stat_name == "I3" ~ "3\nOffspring survival",
                               Stat_name == "I12" ~ "12\nNumber of offspring",
                               Stat_name == "I123" ~ "123\nNumber of offspring\nalive at 60 days",
                               Stat_name == "C1" ~ "Covariance between\nnumber of mates and\nnumber of offspring per mate",
                               Stat_name == "C2" ~ "Covariance between\nnumber of offspring\nand offspring survival",
                               TRUE ~ NA) %>% 
           factor(levels = c("1\nNumber of mates", "2\nNumber of offspring per mate",
                             "Covariance between\nnumber of mates and\nnumber of offspring per mate",
                             "12\nNumber of offspring", "3\nOffspring survival",
                             "Covariance between\nnumber of offspring\nand offspring survival",
                             "123\nNumber of offspring\nalive at 60 days")),
         Percentage = round(Opportunity / I_tot, digits = 1) * 100,
         Sex = ifelse(Sex == "F", "Female", "Male") %>% 
           factor(levels = c("Female", "Male")),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("J. forsmani", "J. praehirsuta")))




########################################
########### Plot the results ###########
########################################
bold_labels <- ifelse(levels(I_vals$Stat_name) %in% c("12\nNumber of offspring", "123\nNumber of offspring\nalive at 60 days"), "bold", "plain")
x_text_colour <- ifelse(levels(I_vals$Stat_name) %in% c("12\nNumber of offspring", "123\nNumber of offspring\nalive at 60 days"), "black", "grey30")

plot_I_vals <- I_vals %>% 
  ggplot() + 
  geom_linerange(aes(x = Stat_name, ymin = CI_2.5, ymax = CI_97.5),
                 lwd = 0.9) + 
  geom_linerange(data = . %>% 
                   filter(Stat_name %in% c("I12",  "I123")),
                 aes(x = Stat_name, ymin = CI_2.5, ymax = CI_97.5),
                 lwd = 1.5) + 
  geom_hline(yintercept = 0, lwd = 1.2, colour = "black") +
  geom_point(aes(x = Stat_name, y = Opportunity, colour = Species),
             size = 4) +
  geom_text(aes(x = Stat_name, y = Inf, label = paste0(Percentage, "%")),
            vjust = 1, size = 5) +
  scale_colour_manual(values = c("J. forsmani" = "#3A9AB2", "J. praehirsuta" = "navy")) +
  ggh4x::facet_grid2(Species ~ Sex, scales = "free") +
  labs(x = "Variance component",
       y = "Value") +
  my_theme +
  theme(strip.text.y = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(face = bold_labels,
                                   colour = x_text_colour,
                                   angle = 90, vjust = 0.5))

plot_I_vals
ggsave("../Figures/tmp/I_values.png", plot = plot_I_vals,
       scale = 3, width = 1400, height = 1500, units = "px")
