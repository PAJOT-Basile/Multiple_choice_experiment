########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "patchwork")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################
########### Useful variables ###########
########################################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

theme_distrib <- my_theme +
  theme(axis.title = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

########################################
############# Import data ##############
########################################
# Survival
survival <- read.table("Associated_data/Raw_metadata/Survival_offspring.tsv",
                       sep = "\t", header = TRUE) %>% 
  as_tibble() %>%
  mutate(Plate_nb = as.numeric(str_remove_all(Plate_ID, "O")))

mothers <- read.table("Associated_data/metadata.tsv",
                      sep = "\t", header = TRUE) %>% 
  as_tibble() %>% 
  filter(Sex == "F", Family_level == "parents",
         !is.na(Species)) %>% 
  unique() %>% 
  select(Label, Species, nb_tot_juv_brood1, nb_tot_juv_brood2, Name_offsprings_brood1, Name_offsprings_brood2) %>% 
  rename(Nbjuvs_1 = nb_tot_juv_brood1, Nbjuvs_2 = nb_tot_juv_brood2, Name_1 = Name_offsprings_brood1, Name_2 = Name_offsprings_brood2) %>%
  pivot_longer(c(c("Nbjuvs_1", "Nbjuvs_2"), c("Name_1", "Name_2")), names_to = c(".value", "Brood"), names_sep = "_") %>% 
  drop_na() %>% 
  filter(Brood == "1") %>% 
  mutate(Min_brood_nb = str_split_fixed(Name, "-", 2)[, 1] %>% str_remove_all("O") %>% as.numeric(),
         Max_brood_nb = str_split_fixed(Name, "-", 2)[, 2] %>% str_remove_all("O") %>% as.numeric(),
         Max_brood_nb = ifelse(is.na(Max_brood_nb), Min_brood_nb, Max_brood_nb))

# Merge the two tables together
merged_table <- survival %>% 
  left_join(mothers, by = join_by("Mother_ID" == "Label"),
            relationship = "many-to-many") %>% 
  filter(Plate_nb >= Min_brood_nb & Plate_nb <= Max_brood_nb) %>%
  unique() %>%
  group_by(Mother_ID, Brood) %>% 
  summarize(Nb_juvs_tot = unique(Nbjuvs),
            Nb_indivs_isolation = sum(Nb_indivs_isolation, na.rm = TRUE),
            Nb_indivs_40d = sum(Nb_indivs_40d, na.rm = TRUE),
            Nb_indivs_60d = sum(Nb_indivs_60d, na.rm = TRUE),
            Species = unique(Species)) %>% 
  mutate(Surv_40 = Nb_indivs_40d / Nb_indivs_isolation,
         Surv_60 = Nb_indivs_60d / Nb_indivs_isolation)



########################################
################ Plots #################
########################################
ff <- merged_table %>%
  filter(Species == "forsmani") %>% 
  ggplot(aes(x = Nb_juvs_tot, y = Surv_60)) +
  geom_point(size = 3, alpha = 0.8, colour = "#3A9AB2") +
  annotate(geom = "text", x = Inf, y = -Inf, label = "(A)",
           size = 8, hjust = 1, vjust = -1.5) +
  xlim(-3, 40) +
  ylim(-0.05, 1.05) +
  labs(x = "Number of offspring per brood\n(W12)",
       y = "Offspring survival\n(W123)") +
  my_theme 

pf <- merged_table %>%
  filter(Species == "praehirsuta") %>% 
  ggplot(aes(x = Nb_juvs_tot, y = Surv_60)) +
  geom_point(size = 3, alpha = 0.8, colour = "navy") +
  annotate(geom = "text", x = Inf, y = -Inf, label = "(B)",
           size = 8, hjust = 1, vjust = -1.5) +
  xlim(-3, 40) +
  ylim(-0.05, 1.05) +
  labs(x = "Number of offspring per brood\n(W12)") +
  my_theme +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())





# Plot the marginal distributions for each fitness component
distrib_w12_ff <- merged_table %>%
  filter(Species == "forsmani") %>% 
  ggplot(aes(x = Nb_juvs_tot)) +
  geom_histogram(fill = "#3A9AB2") +
  xlim(-3, 40) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib

distrib_w12_pf <- merged_table %>%
  filter(Species == "praehirsuta") %>% 
  ggplot(aes(x = Nb_juvs_tot)) +
  geom_histogram(fill = "navy") +
  xlim(-3, 40) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib



distrib_w3_ff <- merged_table %>%
  filter(Species == "forsmani") %>% 
  ggplot(aes(y = Surv_60)) +
  geom_histogram(fill = "#3A9AB2") +
  ylim(-0.05, 1.05) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib

distrib_w3_pf <- merged_table %>%
  filter(Species == "praehirsuta") %>% 
  ggplot(aes(y = Surv_60)) +
  geom_histogram(fill = "navy") +
  ylim(-0.05, 1.05) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib



ann1 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~forsmani)"),
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()
ann2 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~praehirsuta)"),
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()


# Plot everything together
design <- "
A#B#
C#D#
EFGH
"

final_figure <- ann1 + ann2 + 
  distrib_w12_ff + distrib_w12_pf + 
  ff + distrib_w3_ff + pf + distrib_w3_pf +
  plot_layout(design = design,
              widths = c(1, 0.15, 1, 0.15),
              heights = c(0.1, 0.2, 1))

final_figure
ggsave(plot = final_figure, "../Figures/tmp/Distribution_brood_size_survival.png",
       scale = 4, width = 1400, height = 800, units = "px")
