#### Libraries ####
libraries <- c("tidyverse", "ggh4x", "ggfoundry", "ggExtra", "ggpubr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

#### Useful variables ####
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")

set.seed(1234)
#### Useful functions ####
"%not_in%" <- function(x, y) return(!(x %in% y))

str_split_last <- function(string, pattern){
  nb_substrings <- str_count(string, pattern) + 1
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

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
#### Import data ####
# Import metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the dulicates
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  rowwise() %>% 
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all("a|b"),
         Size = str_replace_all(Size, ",", ".") %>% 
           as.numeric(),
         # Recompute the total number of spines and setae per pereopod to homogeneise the values for the parents and offspring
         P7.ep = case_when(
           all(is.na(c(P7.very.large.ep, P7.larg.ep, P7.med.ep, P7.small.ep))) ~ NA,
           TRUE ~ sum(P7.very.large.ep, P7.larg.ep, P7.med.ep, P7.small.ep, na.rm = TRUE)),
         P6.ep.dist = case_when(all(is.na(c(P6.ep.dist.very.large, P6.ep.dist.medium, P6.ep.dist.large, P6.ep.dist.small))) ~ NA,
                                TRUE ~ sum(P6.ep.dist.very.large, P6.ep.dist.medium, P6.ep.dist.large, P6.ep.dist.small, na.rm = TRUE)),
         P6.ep = case_when(
           all(is.na(c(P6.very.large.ep, P6.larg.ep, P6.med.ep, P6.small.ep))) ~ NA,
           TRUE ~ sum(P6.very.large.ep, P6.larg.ep, P6.med.ep, P6.small.ep, na.rm = TRUE)
         ),
         Family_level = factor(Family_level, levels = c("parents", "offspring")),
         Species = factor(Species, levels = c("forsmani", "hybrid", "praehirsuta")),
         # Remove the ">" that are in the data
         across(starts_with(paste0("P", 1:7)), ~ str_remove_all(., ">") %>% trimws() %>% as.numeric())) %>% 
  unique()


# Reproductive success
load("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Sexual_selection/reproductive_success_parents.Rdata")
reproductive_success_parents <- reproductive_success_parents %>% 
  rename(Sex = Parent) %>% 
  mutate(Sex = str_remove_all(Sex, "s")) %>% 
  left_join(metadata %>%
              select(ID_DNA_RAD, Size),
            by = join_by("Parent_ID" == "ID_DNA_RAD")) %>% 
  # Here, we remove the female that had one offspring that was not genotyped
  filter(Parent_ID != "E4961") %>% 
  drop_na() %>% 
  group_by(Species, Sex) %>% 
  rename(W12 = nb_offspring,
         W1 = nb_mates,
         W123 = nb_offspring_surv40,
         W124 = nb_offspring_surv60) %>% 
  mutate(
    # Number of individuals
    p0 = 1/ n(),
    # Standardised fitnesses
    w1 = W1 / mean(W1, na.rm = TRUE),
    w12 = W12 / mean(W12, na.rm = TRUE),
    w123 = W123 / mean(W123, na.rm = TRUE),
    w124 = W124 / mean(W124, na.rm = TRUE),
    # Absolute fitnesses
    W2 = W12 / W1,
    W3 = W123 / W12,
    W4 = W124 / W12,
    # REST
    p1 = p0 * w1,
    w2 = W2 / sum(p1 * W2, na.rm = TRUE),
    p2 = p1 * w2,
    w3 = W3 / sum(p2 * W3, na.rm = TRUE),
    p3 = p2 * w3,
    w4 = W4 / sum(p2 * W4, na.rm = TRUE),
    p4 = p2 * w4,
    across(c(p0, p1, p2,, p3, p4, w1, W2, w2, w12, W12, W3, w3, w123, W123, W4, w4, w124, W124), ~ ifelse(is.na(.), 0, .)))

#### Contributions to total opportunity for selection ####
I_vals <- reproductive_success_parents %>% 
  filter(w1 > 0) %>% 
  group_by(Species, Sex) %>% 
  group_modify(~ bind_rows(coefficients(lm(w12 ~ w1, data = .x)))) %>% 
  select(Species, Sex, w1) %>% 
  rename(Beta_SS = w1) %>% 
  left_join(reproductive_success_parents %>% 
              filter(w1 > 0) %>% 
              group_by(Species, Sex) %>% 
              group_modify(~ bind_rows(confint(lm(w12 ~ w1, data = .x))[2, ] %>% as.data.frame())) %>% 
              rename(Beta_CI = ".") %>% 
              mutate(CI = c("Beta_SS_2.5", "Beta_SS_97.5")) %>% 
              pivot_wider(names_from = "CI", values_from = "Beta_CI"),
            by = c("Species", "Sex")) %>% 
  left_join(reproductive_success_parents %>% 
          group_by(Species, Sex) %>% 
          summarize(Variance = var(W1),
                    I1 = sum(p0 * (w1 - 1)^2)),
          by = c("Species", "Sex")) %>% 
  mutate(smax = Beta_SS * sqrt(I1),
         across(c(contains("I"), contains("Beta"), smax), ~ round(., digits = 2)))

# Compute CIs for the values of I1 and smax
nboots <- 5000
distrib_smax_betass_Is <- lapply(c("forsmani", "praehirsuta"), function(species){
  lapply(c("Male", "Female"), function(sex, sp){
    reproductive_success_parents %>% 
      filter(Sex == sex, Species == sp) %>% 
      bootstrap_smax_betass_Is(., nboots) %>% 
      mutate(Species = sp, Sex = sex)
  }, species)
}) %>% 
  bind_rows()

CI_stats <- distrib_smax_betass_Is %>% 
  replace_na(list("Beta_SS" = 0, "I_s" = 0, "s_max" = 0)) %>% 
  pivot_longer(cols = c(Beta_SS, I_s, s_max), names_to = "Stat_name", values_to = "Stat") %>% 
  group_by(Species, Sex, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat, probs = 0.025),
            CI_97.5 = quantile(Stat, probs = 0.975)) %>% 
  mutate(Stat_name = factor(Stat_name, levels = c("Beta_SS", "I_s", "s_max"),
                            labels = c(expression(Beta[SS]),
                                       expression(I[S]),
                                       expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')")))

calculated_values <- I_vals %>%
  pivot_longer(c(Beta_SS, I1, smax), names_to = "Stat_name", values_to = "Stat") %>% 
  mutate(Stat_name = case_when(Stat_name == "I1" ~ "I_s",
                               Stat_name == "smax" ~ "s_max",
                               TRUE ~ Stat_name) %>% 
           factor(levels = c("Beta_SS", "I_s", "s_max"),
                  labels = c(expression(Beta[SS]),
                             expression(I[S]),
                             expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')"))) 

distrib_params <- distrib_smax_betass_Is %>%
  replace_na(list("Beta_SS" = 0, "I_s" = 0, "s_max" = 0)) %>% 
  pivot_longer(cols = c(Beta_SS, I_s, s_max), names_to = "Stat_name", values_to = "Stat") %>% 
  mutate(Stat_name = factor(Stat_name, levels = c("Beta_SS", "I_s", "s_max"),
                            labels = c(expression(Beta[SS]),
                                       expression(I[S]),
                                       expression("s'"[max] == Beta[SS]~sqrt(I[S])))),
         Species = factor(Species, levels = c("forsmani", "praehirsuta"),
                          labels = c("italic('J. forsmani')", "italic('J. praehirsuta')"))) %>% 
  ggplot() +
  geom_rect(data = CI_stats, aes(xmin = CI_2.5, xmax = CI_97.5, ymin = -Inf, ymax = Inf, fill = Sex), alpha = 0.2) +
  geom_histogram(aes(x = Stat, fill = Sex), bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(data = calculated_values, aes(xintercept = Stat, colour = Sex), lwd = 1, lty = 2) +
  scale_fill_manual(values = colours_sexes) +
  scale_colour_manual(values = colours_sexes) +
  labs(x = "Statistic value",
       y = "Count") +
  facet_grid2(Species ~ Stat_name, labeller = as_labeller(label_parsed), scales = "free") +
  my_theme

distrib_params
ggsave("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Sexual_selection/Distribution_betaSS_Is_smax.png",
       distrib_params, width = 1200, height = 750, units = "px", dpi = 500, scale = 8)
