########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "ggh4x", "ggfoundry", "ggExtra", "ggpubr", "patchwork")
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

set.seed(1234)
########################################
########### Useful functions ###########
########################################
"%not_in%" <- function(x, y) return(!(x %in% y))

compute_variances_in_bootstrap <- function(df){
  # Sample randomly columns
  lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  df_samp <- df[lines_to_sample, ]
  
  # Calculate the parameters
  df_samp$w1 <- df_samp$W1 / mean(df_samp$W1, na.rm = TRUE)
  df_samp$w12 <- df_samp$W12 / mean(df_samp$W12, na.rm = TRUE)
  df_samp$w123 <- df_samp$W123 / mean(df_samp$W123, na.rm = TRUE)
  df_samp$w1234 <- df_samp$W1234 / mean(df_samp$W1234, na.rm = TRUE)
  df_samp$w1235 <- df_samp$W1235 / mean(df_samp$W1235, na.rm = TRUE)
  
  df_samp$p1 <- df_samp$p0 * df_samp$w1
  df_samp$w2 <- df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE)
  
  df_samp$p2 <- df_samp$p1 * df_samp$w2
  df_samp$w3 <- df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE)
  
  df_samp$p3 <- df_samp$p2 * df_samp$w3
  df_samp$w4 <- df_samp$W4 / sum(df_samp$p3 * df_samp$W4, na.rm = TRUE)
  df_samp$w5 <- df_samp$W5 / sum(df_samp$p3 * df_samp$W5, na.rm = TRUE)
  df_samp <- df_samp %>% 
    mutate(across(c(w1, w12, w123, w1234, w1235, p1, w2, p2, w3, p3, w4, w5), ~ ifelse(is.na(.), 0, .)))

  # Calculate the statistics
  I1 <- sum(df_samp$p0 * (df_samp$w1 - 1)^2, na.rm = TRUE)
  I2 <- sum(df_samp$p1 * (df_samp$w2 - 1)^2, na.rm = TRUE)
  I12 <- sum(df_samp$p0 * (df_samp$w12 - 1)^2, na.rm = TRUE)
  
  COI1.2 <- sum(df_samp$p0 * (df_samp$w1 - 1) * (df_samp$w2 - 1), na.rm = TRUE)
  COI1.2_1 <- sum(df_samp$p1 * (df_samp$w1 - 1) * (df_samp$w2 - 1), na.rm = TRUE)
  COI12.2_1 <- sum(df_samp$p1 * (df_samp$w12 - 1) * (df_samp$w2 - 1), na.rm = TRUE)
  COI12.2 <- sum(df_samp$p0 * (df_samp$w12 - 1) * (df_samp$w2 - 1), na.rm = TRUE)
  
  I3 <- sum(df_samp$p2 * (df_samp$w3 - 1)^2, na.rm = TRUE)
  I123 <- sum(df_samp$p0 * (df_samp$w123 - 1)^2, na.rm = TRUE)
  
  COI1.3 <- sum(df_samp$p0 * (df_samp$w1 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  COI2.3 <- sum(df_samp$p1 * (df_samp$w2 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  COI12.3 <- sum(df_samp$p0 * (df_samp$w12 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  COI12.3_2 <- sum(df_samp$p2 * (df_samp$w12 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  COI123.3_2 <- sum(df_samp$p2 * (df_samp$w123 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  COI123.3 <- sum(df_samp$p0 * (df_samp$w123 - 1) * (df_samp$w3 - 1), na.rm = TRUE)
  
  I4 <- sum(df_samp$p3 * (df_samp$w4 - 1)^2, na.rm = TRUE)
  I1234 <- sum(df_samp$p0 * (df_samp$w1234 - 1)^2, na.rm = TRUE)
  
  COI1.4 <- sum(df_samp$p0 * (df_samp$w1 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI2.4 <- sum(df_samp$p1 * (df_samp$w2 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI3.4 <- sum(df_samp$p2 * (df_samp$w3 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI123.4 <- sum(df_samp$p0 * (df_samp$w123 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI123.4_3 <- sum(df_samp$p3 * (df_samp$w123 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI1234.4_3 <- sum(df_samp$p3 * (df_samp$w1234 - 1) * (df_samp$w4 - 1), na.rm = TRUE)
  COI1234.4 <- sum(df_samp$p0 * (df_samp$w1234 - 1) * (df_samp$w4 - 1), na.rm = TRUE)

  I5 <- sum(df_samp$p3 * (df_samp$w5 - 1)^2, na.rm = TRUE)
  
  COI1.5 <- sum(df_samp$p0 * (df_samp$w1 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI2.5 <- sum(df_samp$p1 * (df_samp$w2 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI3.5 <- sum(df_samp$p2 * (df_samp$w3 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI123.5 <- sum(df_samp$p0 * (df_samp$w123 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI123.5_3 <- sum(df_samp$p3 * (df_samp$w123 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI1235.5_3 <- sum(df_samp$p3 * (df_samp$w1235 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  COI1235.5 <- sum(df_samp$p0 * (df_samp$w1235 - 1) * (df_samp$w5 - 1), na.rm = TRUE)
  
  I_tot <- sum(df_samp$p0 * (df_samp$w1235 - 1)^2, na.rm = TRUE)
  
  to_ret <- c("I1" = I1, "I2" = I2, "I12" = I12,
              "COI1.2" = COI1.2, "COI1.2_1" = COI1.2_1, "COI12.2_1" = COI12.2_1, "COI12.2" = COI12.2,
              "I3" = I3, "I123" = I123,
              "COI1.3" = COI1.3, "COI2.3" = COI2.3, "COI12.3" = COI12.3, "COI12.3_2" = COI12.3_2, "COI123.3_2" = COI123.3_2, "COI123.3" = COI123.3,
              "I4" = I4, "I1234" = I1234,
              "COI1.4" = COI1.4, "COI2.4" = COI2.4, "COI3.4" = COI3.4, "COI123.4" = COI123.4, "COI123.4_3" = COI123.4_3, "COI1234.4_3" = COI1234.4_3, "COI1234.4" = COI1234.4,
              "I5" = I5,
              "COI1.5" = COI1.5, "COI2.5" = COI2.5, "COI3.5" = COI3.5, "COI123.5" = COI123.5, "COI123.5_3" = COI123.5_3, "COI1235.5_3" = COI1235.5_3, "COI1235.5" = COI1235.5,
              "I_tot" = I_tot)
  return(to_ret)
}

bootstrap_variances <- function(df, nb_boots){
  replicate(nb_boots, 
            expr = {compute_variances_in_bootstrap(df)},
            simplify = "vector") %>% 
    t() %>% 
    as_tibble() %>% 
    return()
}


########################################
############# Import data ##############
########################################
# Import metadata
metadata <- read.table("Associated_data//metadata.tsv",
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
load("Associated_data/Rdata/reproductive_success_parents.Rdata")
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
  rename(W1 = nb_mates,
         W12 = max_nb_offspring,
         W123 = nb_offspring,
         W1234 = nb_offspring_surv40,
         W1235 = nb_offspring_surv60) %>% 
  mutate(
    # Number of individuals
    p0 = 1/ n(),
    # Standardised fitnesses
    w1 = W1 / mean(W1, na.rm = TRUE),
    w12 = W12 / mean(W12, na.rm = TRUE),
    w123 = W123 / mean(W123, na.rm = TRUE),
    w1234 = W1234 / mean(W1234, na.rm = TRUE),
    w1235 = W1235 / mean(W1235, na.rm = TRUE),
    # Absolute fitnesses
    W2 = W12 / W1,
    W3 = W123 / W12,
    W4 = W1234 / W123,
    W5 = W1235 / W123,
    # REST
    p1 = p0 * w1,
    w2 = W2 / sum(p1 * W2, na.rm = TRUE),
    p2 = p1 * w2,
    w3 = W3 / sum(p2 * W3, na.rm = TRUE),
    p3 = p2 * w3,
    w4 = W4 / sum(p3 * W4, na.rm = TRUE),
    w5 = W5 / sum(p3 * W5, na.rm = TRUE),
    across(c(w1, w12, w123, w1234, w1235, p1, w2, p2, w3, p3, w4, w5), ~ ifelse(is.na(.), 0, .)))



########################################
## Compute opportunity for selection ###
########################################
I_vals <- reproductive_success_parents %>% 
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
                        
                        I4 = sum(p3 * (w4 - 1)^2, na.rm = TRUE),
                        I1234 = sum(p0 * (w1234 - 1)^2, na.rm = TRUE),
                        
                        COI1.4 = sum(p0 * (w1 - 1) * (w4 - 1), na.rm = TRUE),
                        COI2.4 = sum(p1 * (w2 - 1) * (w4 - 1), na.rm = TRUE),
                        COI3.4 = sum(p2 * (w3 - 1) * (w4 - 1), na.rm = TRUE),
                        COI123.4 = sum(p0 * (w123 - 1) * (w4 - 1), na.rm = TRUE),
                        COI123.4_3 = sum(p3 * (w123 - 1) * (w4 - 1), na.rm = TRUE),
                        COI1234.4_3 = sum(p3 * (w1234 - 1) * (w4 - 1), na.rm = TRUE),
                        COI1234.4 = sum(p0 * (w1234 - 1) * (w4 - 1), na.rm = TRUE),
                        
                        I5 = sum(p3 * (w5 - 1)^2, na.rm = TRUE),
                        
                        COI1.5 = sum(p0 * (w1 - 1) * (w5 - 1), na.rm = TRUE),
                        COI2.5 = sum(p1 * (w2 - 1) * (w5 - 1), na.rm = TRUE),
                        COI3.5 = sum(p2 * (w3 - 1) * (w5 - 1), na.rm = TRUE),
                        COI123.5 = sum(p0 * (w123 - 1) * (w5 - 1), na.rm = TRUE),
                        COI123.5_3 = sum(p3 * (w123 - 1) * (w5 - 1), na.rm = TRUE),
                        COI1235.5_3 = sum(p3 * (w1235 - 1) * (w5 - 1), na.rm = TRUE),
                        COI1235.5 = sum(p0 * (w1235 - 1) * (w5 - 1), na.rm = TRUE),
                        
                        I_tot = sum(p0 * (w1235 - 1)^2, na.rm = TRUE)) %>% 
  mutate(C1 = COI1.2 + COI1.2_1 + (COI12.2_1 - COI12.2),
         C2 = COI12.3 + COI12.3_2 + (COI123.3_2 - COI123.3),
         C3 = COI123.4 + COI123.4_3 + (COI1234.4_3 - COI1234.4),
         C4 = COI123.5 + COI123.5_3 + (COI1235.5_3 - COI1235.5),
         across(c(contains("I", ignore.case = FALSE), contains("C")),
                ~ (. / I_tot) * 100,
                .names = "Percentage_{col}"),
         across(c(contains("I", ignore.case = FALSE), contains("C")), ~ round(., digits = 2)))

########################################
######### Confidence intervals #########
########################################
# nboots <- 50000
# distrib_variances <- lapply(c("forsmani", "praehirsuta"), function(species){
#   lapply(c("Male", "Female"), function(sex, sp){
#     reproductive_success_parents %>% 
#       filter(Sex == sex, Species == sp) %>% 
#       bootstrap_variances(., nboots) %>% 
#       mutate(Species = sp, Sex = sex)
#   }, species)
# }) %>% 
#   bind_rows() %>% 
#   mutate(C1 = COI1.2 + COI1.2_1 + (COI12.2_1 - COI12.2),
#          C2 = COI12.3 + COI12.3_2 + (COI123.3_2 - COI123.3),
#          C3 = COI123.4 + COI123.4_3 + (COI1234.4_3 - COI1234.4),
#          C4 = COI123.5 + COI123.5_3 + (COI1235.5_3 - COI1235.5))
# 
# save(distrib_variances, file = "Associated_data/Rdata/distribution_variances_50000_boots.rda")
load("Associated_data/Rdata/distribution_variances_50000_boots.rda")


CI_stats <- distrib_variances %>% 
  pivot_longer(cols = !c(Sex, Species), names_to = "Stat_name", values_to = "Stat") %>% 
  group_by(Species, Sex, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat, probs = 0.025),
            CI_97.5 = quantile(Stat, probs = 0.975))

I_vals %>% 
  pivot_longer(!c(Species, Sex), names_to = "Stat_name", values_to = "Stat") %>% 
  left_join(CI_stats, by = c("Species", "Sex", "Stat_name")) %>% 
  rowwise() %>% 
  mutate(CI = paste0("[", round(CI_2.5, digits = 2), "; ", round(CI_97.5, digits = 2), "]")) %>% 
  select(-c(CI_2.5, CI_97.5)) %>% 
  filter(Stat_name == "COI1.2") 

calculated_values <- I_vals %>% 
  pivot_longer(cols = !c(Sex, Species), names_to = "Stat_name", values_to = "Stat")

########################################
########## Bateman's gradient ##########
########################################
# Plot bateman's gradient for both sexes of both species
fm <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Male") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "#3A9AB2") +
  geom_smooth(data = reproductive_success_parents %>% 
                filter(Species == "forsmani", Sex == "Male", w1 > 0),
              colour = "#3A9AB2", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(B)",
           size = 10, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  my_theme +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


ff <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Female") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "#3A9AB2") +
  geom_smooth(data = reproductive_success_parents %>% 
                filter(Species == "forsmani", Sex == "Female", w1 > 0),
              colour = "#3A9AB2", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(A)",
           size = 10, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  my_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

pm <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
                Sex == "Male") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "navy") +
  geom_smooth(data = reproductive_success_parents %>% 
                filter(Species == "praehirsuta", Sex == "Male", w1 > 0),
              colour = "navy", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(D)",
           size = 10, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pf <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
         Sex == "Female") %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_point(size = 3, alpha = 0.8, colour = "navy") +
  geom_smooth(data = reproductive_success_parents %>% 
                filter(Species == "praehirsuta", Sex == "Female", w1 > 0),
              colour = "navy", method = "lm", se = FALSE) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "(C)",
           size = 10, hjust = -.5, vjust = 1.5) +
  xlim(-0.1, 3) +
  ylim(-0.107, 4) +
  my_theme +
  theme(axis.title = element_blank())


# Plot the marginal distributions for each fitness component
distrib_w1_ff <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Female") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "#3A9AB2") +
  xlim(-0.1, 3) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w12_ff <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Female") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "#3A9AB2") +
  ylim(-2, 75) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w1_fm <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Male") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "#3A9AB2") +
  xlim(-0.1, 3) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w12_fm <- reproductive_success_parents %>%
  filter(Species == "forsmani",
         Sex == "Male") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "#3A9AB2") +
  ylim(-2, 75) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w1_pf <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
         Sex == "Female") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "navy") +
  xlim(-0.1, 3) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w12_pf <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
         Sex == "Female") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "navy") +
  ylim(-2, 75) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w1_pm <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
         Sex == "Male") %>% 
  ggplot(aes(x = W1)) +
  geom_histogram(fill = "navy") +
  xlim(-0.1, 3) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

distrib_w12_pm <- reproductive_success_parents %>%
  filter(Species == "praehirsuta",
         Sex == "Male") %>% 
  ggplot(aes(y = W12)) +
  geom_histogram(fill = "navy") +
  ylim(-2, 75) +
  my_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank())

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

ann5 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Relative~reproductive~success~(w12))"),
            parse = TRUE, size = 6, hjust = 0.5, angle = 90) +
  theme_void()

ann6 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Relative~mating~success~(w1))"),
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()


# Plot everything together
design <- "
#A#B##
#C#D##
EFGHIJ
EK#L##
EMNOPQ
#RRR##
"

final_figure <- ann1 + ann2 + 
  distrib_w1_ff + distrib_w1_fm + 
  ann5 + ff + distrib_w12_ff + fm + distrib_w12_fm + ann3 + 
  distrib_w1_pf + distrib_w1_pm +
  pf + distrib_w12_pf + pm + distrib_w12_pm + ann4 + 
  ann6 +
  plot_layout(design = design,
              widths = c(0.1, 1, 0.15, 1, 0.15, 0.1),
              heights = c(0.1, 0.2, 1, 0.2, 1, 0.1))

final_figure
ggsave(plot = final_figure, "../Figures/Distribution_reproductive_success.png", scale = 4, width = 1200, height = 1100, units = "px")
