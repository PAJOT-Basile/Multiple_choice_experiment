# Compute confidence intervals for different sexual selection parameters
# and plot the distributions of their bootstraps:
#   1. Bateman's gradient
#   2. Opportunity for sexual selection
#   3. Jones's index s'max

# Finally, plot Bateman's gradient for the different sexes of the two species

########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "ggh4x", "patchwork")
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

colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")

set.seed(1234)

nboots <- 50000

########################################
########### Useful functions ###########
########################################
source("98_Transversal_scripts/Transversal_functions.r")

compute_smax <- function(beta_SS, I_s){
  return(beta_SS * sqrt(I_s))
}

compute_Is <- function(df){
  return(sum(df$p0 * (df$w1 - 1)^2))
}

compute_BetaSS <- function(df){
  return(lm(w12 ~ w1, data = df[which(df$w1 > 0), ])$coef[2])
}

compute_smax_betass_Is_in_bootstrap <- function(df){
  lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  df_samp <- df[lines_to_sample, ]
  df_samp$w1 <- df_samp$W1 / mean(df_samp$W1, na.rm = TRUE)
  df_samp$w12 <- df_samp$W12 / mean(df_samp$W12, na.rm = TRUE)
  beta_SS <- compute_BetaSS(df_samp)
  I_s <- compute_Is(df_samp)
  s_max <- compute_smax(beta_SS, I_s)
  
  return(c("Beta_SS" = beta_SS, "I_s" = I_s, "s_max" = s_max))
}

bootstrap_smax_betass_Is <- function(df, nb_boots){
  replicate(nb_boots, 
            expr = {compute_smax_betass_Is_in_bootstrap(df)},
            simplify = "vector") %>% 
    t() %>% 
    as_tibble() %>% 
    rename(Beta_SS = Beta_SS.w1,
           s_max = s_max.w1) %>% 
    return()
}


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
  group_by(Species, Sex) %>% 
  rename(W1 = nb_partners,
         W12 = nb_offspring,
         W123 = nb_offspring_surv40,
         W124 = nb_offspring_surv60) %>% 
  mutate(
    # Number of individuals
    p0 = 1/ n(),
    # Standardised fitnesses
    w1 = W1 / mean(W1, na.rm = TRUE),
    w12 = W12 / mean(W12, na.rm = TRUE),
    across(c(p0, w1, w12, W12), ~ ifelse(is.na(.), 0, .)))

interesting_stats <- fitness_components %>% 
  select(Species, Sex, W1, W12) %>% 
  pivot_longer(starts_with("W"), names_to = "Episode", values_to = "Fitness") %>% 
  group_by(Species, Sex, Episode) %>% 
  summarize(min = round(min(Fitness, na.rm = TRUE), digits = 2),
            mean = round(mean(Fitness, na.rm = TRUE), digits = 2),
            max = round(max(Fitness, na.rm = TRUE), digits = 2),
            var = round(var(Fitness, na.rm = TRUE), digits = 2))


#########################################
## Compute sexual selection statistics ##
#########################################
Sexual_selection_stats <- fitness_components %>% 
  filter(w1 > 0) %>% 
  group_by(Species, Sex) %>% 
  group_modify(~ bind_rows(coefficients(lm(w12 ~ w1, data = .x)))) %>% 
  select(Species, Sex, w1) %>% 
  rename(Beta_SS = w1) %>% 
  left_join(fitness_components %>% 
              filter(w1 > 0) %>% 
              group_by(Species, Sex) %>% 
              group_modify(~ bind_rows(summary(lm(w12 ~ w1, data = .x))$coefficients %>%
                                         as.data.frame() %>% 
                                         rownames_to_column("param") %>% 
                                         filter(param == "w1") %>% 
                                         select("Pr(>|t|)"))) %>% 
              rename(pval_beta = "Pr(>|t|)"),
            by = c("Species", "Sex")) %>% 
  left_join(fitness_components %>% 
          group_by(Species, Sex) %>% 
          summarize(I1 = sum(p0 * (w1 - 1)^2)),
          by = c("Species", "Sex")) %>% 
  mutate(smax = Beta_SS * sqrt(I1),
         across(c(contains("I"), contains("Beta"), smax), ~ round(., digits = 2)))

#########################################
###### Compute confidence intervals #####
#########################################
# distrib_smax_betass_Is <- lapply(c("forsmani", "praehirsuta"), function(species){
#   lapply(c("M", "F"), function(sex, sp){
#     fitness_components %>%
#       filter(Sex == sex, Species == sp) %>%
#       bootstrap_smax_betass_Is(., nboots) %>%
#       mutate(Species = sp, Sex = sex)
#   }, species)
# }) %>%
#   bind_rows()
# 
# save(distrib_smax_betass_Is, file = "Associated_data/Rdata/distribution_smax_betass_Is_50000_bootstraps.rda")
load("Associated_data/Rdata/distribution_smax_betass_Is_50000_bootstraps.rda")

CI_stats <- distrib_smax_betass_Is %>% 
  replace_na(list("Beta_SS" = 0, "I_s" = 0, "s_max" = 0)) %>% 
  pivot_longer(cols = c(Beta_SS, I_s, s_max), names_to = "Stat_name", values_to = "Stat") %>% 
  group_by(Species, Sex, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat, probs = 0.025, na.rm = FALSE),
            CI_97.5 = quantile(Stat, probs = 0.975, na.rm = FALSE)) %>% 
  # For J. forsmani, there is no variability on w1, so the calculated values in the distribution
  # of the betaSS are all NA. The confidence interval is thus NA.
  mutate(across(c(CI_2.5, CI_97.5), ~ case_when(Sex == "F" & Species == "forsmani" & Stat_name %in% c("Beta_SS", "s_max") ~ NA,
                                                TRUE ~ .))) %>% 
  mutate(Stat_name = factor(Stat_name, levels = c("Beta_SS", "I_s", "s_max"),
                            labels = c(expression(Beta[SS]),
                                       expression(I[S]),
                                       expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')")))

calculated_values <- Sexual_selection_stats %>%
  pivot_longer(c(Beta_SS, I1, smax), names_to = "Stat_name", values_to = "Stat") %>% 
  mutate(Stat_name = case_when(Stat_name == "I1" ~ "I_s",
                               Stat_name == "smax" ~ "s_max",
                               TRUE ~ Stat_name) %>% 
           factor(levels = c("Beta_SS", "I_s", "s_max"),
                  labels = c(expression(Beta[SS]),
                             expression(I[S]),
                             expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')"))) %>% 
  left_join(CI_stats, by = c("Stat_name", "Sex", "Species")) %>% 
  mutate(Recap = case_when(is.na(Stat) ~ NA,
                           TRUE ~ paste0(Stat, " [", round(CI_2.5, digits = 2), "; ", round(CI_97.5, digits = 2), "]")),
         pval_beta = case_when(Stat_name == expression(Beta[SS]) ~ pval_beta,
                               TRUE ~ NA))

distrib_params <- distrib_smax_betass_Is %>%
  # replace_na(list("Beta_SS" = 0, "I_s" = 0, "s_max" = 0)) %>% 
  pivot_longer(cols = c(Beta_SS, I_s, s_max), names_to = "Stat_name", values_to = "Stat") %>% 
  mutate(Stat_name = factor(Stat_name, levels = c("Beta_SS", "I_s", "s_max"),
                            labels = c(expression(Beta[SS]),
                                       expression(I[S]),
                                       expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')")),
         Sex = ifelse(Sex == "F", "Female", "Male")) %>% 
  ggplot() +
  geom_rect(data = CI_stats %>% 
              mutate(Sex = ifelse(Sex == "F", "Female", "Male")),
            aes(xmin = CI_2.5, xmax = CI_97.5, ymin = -Inf, ymax = Inf, fill = Sex), alpha = 0.2) +
  geom_histogram(aes(x = Stat, fill = Sex), bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(data = calculated_values, aes(xintercept = Stat, colour = Sex), lwd = 1, lty = 2) +
  scale_fill_manual(values = colours_sexes) +
  scale_colour_manual(values = colours_sexes) +
  labs(x = "Statistic value",
       y = "Count") +
  facet_grid2(Species ~ Stat_name, labeller = as_labeller(label_parsed), scales = "free") +
  my_theme

distrib_params
ggsave("../Figures/Suppfigs/Distribution_betaSS_Is_smax.png",
       distrib_params, width = 1200, height = 750, units = "px", dpi = 500, scale = 8)


########################################
########## Bateman's gradient ##########
########################################
# Plot bateman's gradient for both sexes of both species
ff <- fitness_components %>%
  filter(Species == "forsmani",
         Sex == "F") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "#3A9AB2") +
  geom_smooth(data = fitness_components %>% 
                filter(Species == "forsmani", Sex == "Female", w1 > 0),
              colour = "#3A9AB2", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(A)",
           size = 8, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  labs(y = "Relative fecundity\n(w12)") +
  my_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

fm <- fitness_components %>%
  filter(Species == "forsmani",
         Sex == "M") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "#3A9AB2") +
  geom_smooth(data = fitness_components %>% 
                filter(Species == "forsmani", Sex == "Male", w1 > 0),
              colour = "#3A9AB2", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(B)",
           size = 8, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  my_theme +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


pf <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "F") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "navy") +
  geom_smooth(data = fitness_components %>% 
                filter(Species == "praehirsuta", Sex == "Female", w1 > 0),
              colour = "navy", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(C)",
           size = 8, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  labs(x = "Relative mating success\n(w1)",
       y = "Relative fecundity\n(w12)") +
  my_theme

pm <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "M") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "navy") +
  geom_smooth(data = fitness_components %>% 
                filter(Species == "praehirsuta", Sex == "Male", w1 > 0),
              colour = "navy", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(D)",
           size = 8, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  labs(x = "Relative mating success\n(w1)") +
  my_theme +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




# Plot the marginal distributions for each fitness component
distrib_w1_ff <- fitness_components %>%
  filter(Species == "forsmani",
         Sex == "F") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "#3A9AB2") +
  xlim(-0.1, 3) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib

distrib_w1_fm <- fitness_components%>%
  filter(Species == "forsmani",
         Sex == "M") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "#3A9AB2") +
  xlim(-0.1, 3) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib

distrib_w1_pf <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "F") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "navy") +
  xlim(-0.1, 3) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib


distrib_w1_pm <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "M") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "navy") +
  xlim(-0.1, 3) +
  scale_y_continuous(n.breaks = 3) +
  theme_distrib

distrib_w12_ff <- fitness_components %>%
  filter(Species == "forsmani",
         Sex == "F") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "#3A9AB2") +
  ylim(-2, 75) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib

distrib_w12_fm <- fitness_components %>%
  filter(Species == "forsmani",
         Sex == "M") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "#3A9AB2") +
  ylim(-2, 75) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib

distrib_w12_pf <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "F") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "navy") +
  ylim(-2, 75) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib


distrib_w12_pm <- fitness_components %>%
  filter(Species == "praehirsuta",
         Sex == "M") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "navy") +
  ylim(-2, 75) +
  scale_x_continuous(n.breaks = 3) +
  theme_distrib



# Prepare the annotations
ann1 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Females)"),
            parse = TRUE, size = 6, hjust = 0.1) +
  theme_void()
ann2 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Males)"),
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()

ann3 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~forsmani)"),
            parse = TRUE, size = 6, hjust = 0.5, angle = 90) +
  theme_void()
ann4 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~praehirsuta)"),
            parse = TRUE, size = 6, hjust = 0.5, angle = 90) +
  theme_void()


# Plot everything together
design <- "
#A#B#
#C#D#
EFGHI
#J#K#
LMNOP
"

final_figure <- ann1 + ann2 + 
  distrib_w1_ff + distrib_w1_fm + 
  ann3 + ff + distrib_w12_ff + fm + distrib_w12_fm +
  distrib_w1_pf + distrib_w1_pm +
  ann4 + pf + distrib_w12_pf + pm + distrib_w12_pm +
  plot_layout(design = design,
              widths = c(0.1, 1, 0.15, 1, 0.15),
              heights = c(0.1, 0.2, 1, 0.2, 1))

final_figure
ggsave(plot = final_figure, "../Figures/Distribution_reproductive_success.png", scale = 4, width = 1200, height = 1100, units = "px")