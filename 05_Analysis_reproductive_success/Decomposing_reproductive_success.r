#### Libraries ####
libraries <- c("tidyverse", "ggh4x", "ggfoundry", "ggExtra", "ggpubr")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

#### Useful variables ####
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("forsmani" = "#3A9AB2", "praehirsuta" = "navy")

set.seed(1234)
#### Useful functions ####
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
load("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Rdata/reproductive_success_parents.Rdata")
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




#### Contributions to total opportunity for selection ####
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

# Compute confidence intervals for the statistics
nboots <- 50000
distrib_variances <- lapply(c("forsmani", "praehirsuta"), function(species){
  lapply(c("Male", "Female"), function(sex, sp){
    reproductive_success_parents %>% 
      filter(Sex == sex, Species == sp) %>% 
      bootstrap_variances(., nboots) %>% 
      mutate(Species = sp, Sex = sex)
  }, species)
}) %>% 
  bind_rows() %>% 
  mutate(C1 = COI1.2 + COI1.2_1 + (COI12.2_1 - COI12.2),
         C2 = COI12.3 + COI12.3_2 + (COI123.3_2 - COI123.3),
         C3 = COI123.4 + COI123.4_3 + (COI1234.4_3 - COI1234.4),
         C4 = COI123.5 + COI123.5_3 + (COI1235.5_3 - COI1235.5))

distrib_variances %>% ggplot(aes(x = COI1.3)) + geom_histogram(bins = 100) + facet_grid2(Species ~ Sex)

CI_stats <- distrib_variances %>% 
  # select(Sex, Species, I1, I2, C1, I12, I4, C2, I_tot) %>% 
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
  # select(Sex, Species, I1, I2, C1, I12, I4, C2, I_tot) %>% 
  pivot_longer(cols = !c(Sex, Species), names_to = "Stat_name", values_to = "Stat")

######### Look at distributions  #########
# W1
reproductive_success_parents %>% 
  # filter(Parent_ID %in% c("E4911", "E4929", "E4940")) %>% 
  # mutate(W1 = as.factor(W1)) %>% 
  
  ggplot(aes(x = W1, y = W5, colour = W1)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(aes(label = Parent_ID)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_gradientn(colors=c("yellow", "blue", "red")) +
  # ggplot() +
  # geom_histogram(aes(x = W5, fill = Species, colour = Species), bins = 20, alpha = 0.5) +
  # scale_fill_manual(values = colours_species) +
  # scale_colour_manual(values = colours_species) +
  facet_grid(Species ~ Sex, scales = "free") +
  my_theme



I_vals %>%
  select(contains("Percentage")) %>%
  pivot_longer(contains("Percentage"), names_to = "Stat_name", values_to = "Stat") %>%
  mutate(Stat_name = str_remove_all(Stat_name, "Percentage_") %>% 
           factor(levels = c("I1", "I2", "COI1.2", "COI1.2_1", "change_cov_nb_offs_fertility", "I12", "I4", "COI12.4", "COI12.4_2", "change_cov_nb_offs_surv", "I_tot"))) %>%
  ggplot(aes(x = Stat_name, y = Stat)) +
  geom_col() +
  facet_grid2(Sex ~ Species) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))


##### Batemans' gradient #####
# Plot
reproductive_success_parents %>% 
  # filter(Parent_ID %not_in% c("E4895", "E4925")) %>% 
  ggplot(aes(x = w1, y = w12, colour = Species)) +
  geom_jitter(size = 2, width = 0.04, height = 0) +
  geom_smooth(data = reproductive_success_parents %>% 
                filter(w1 != 0),
              method = "lm", se = FALSE) +
  scale_colour_manual(values = c("forsmani" = "#3A9AB2", "praehirsuta" = "navy")) +
  facet_grid2(Species ~ Sex,
              labeller = as_labeller(c(
                "forsmani" = "J. forsmani", "praehirsuta" = "J. praehirsuta",
                "Female" = "Female",
                "Male" = "Male"))) +
  labs(x = "Relative mating success",
       y = "Relative reproductive success") +
  my_theme +
  theme(strip.text.y = element_text(face = "italic"))




# Represent I1
reproductive_success_parents %>%
  left_join(I_vals, by = c("Species", "Sex")) %>% 
  ggplot(aes(x = w1)) +
  geom_histogram(fill = "white", colour = "black") +
  geom_segment(aes(x = mean(w1) - var(w1)/2,
                   xend = mean(w1) + var(w1)/2,
                   y = 5),
               lwd = 0.05,
               arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_segment(aes(xend = mean(w1) - var(w1)/2,
                   x = mean(w1) + var(w1)/2,
                   y = 5),
               lwd = 0.05,
               arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_text(aes(x = mean(w1) + .01, y = 6, label = paste("I1 =", I1)), size = 3.5) +
  facet_grid2(Species ~ Sex, scales = "free") +
  my_theme +
  labs(x = "Relative number of mates",
       y = "Count")

# Correlation size and relative number of mates
reproductive_success_parents %>% 
  ggplot(aes(x = Size, y = w1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  my_theme + 
  labs() +
  facet_grid2(Species ~ Sex, scales = "free")

# Correlation relative number of mates and phenotypic traits
reproductive_success_parents %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species, P7.ep, P7.ep.mid, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
            by = join_by("Parent_ID" == "ID_DNA_RAD", "Species")) %>% 
  pivot_longer(cols = c(Size, P7.ep, P7.ep.mid, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
               values_to = "count", names_to = "Trait") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, y = w1)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_count() +
  facet_grid2(Species ~ Trait, scales = "free") +
  my_theme +
  labs(x = "Number of curved setae/spines per pereiopod",
       y = "Standardised number of mates")


analyse_directional_selection <- reproductive_success_parents %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species, P7.ep, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
            by = join_by("Parent_ID" == "ID_DNA_RAD", "Species")) %>% 
  pivot_longer(cols = c(Size, P7.ep, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
               values_to = "z",
               names_to = "Trait") %>% 
  drop_na() %>% 
  group_by(Species, Sex, Trait) %>% 
  summarize(z1_bar = sum(p1 * z, na.rm = TRUE),
            z_bar = sum(p0 * z, na.rm = TRUE),
            z2_bar = sum(p2 * z, na.rm = TRUE),
            S1 = sum(p1 * z, na.rm = TRUE) - sum(p0 * z, na.rm = TRUE),
            i1 = S1 / sd(z, na.rm = TRUE),
            beta_1 = S1 / var(z, na.rm = TRUE),
            S2 = sum(p2 * z, na.rm = TRUE) - sum(p1 * z, na.rm = TRUE),
            i2 = S2 / sd(z, na.rm = TRUE),
            beta_2 = S2 / var(z, na.rm = TRUE))



reproductive_success_parents %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species, P7.ep, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
            by = join_by("Parent_ID" == "ID_DNA_RAD", "Species")) %>% 
  pivot_longer(cols = c(Size, P7.ep, P6.ep, P6.ep.dist, P5.ep, P5.curv.setae, P4.small.ep, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae),
               values_to = "count", names_to = "Trait") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, y = w1)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_count() +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female"),
             aes(xintercept = z1_bar)) +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female"),
             aes(xintercept = z_bar), colour = "blue") +
  facet_grid2(Species ~ Trait, scales = "free") +
  my_theme +
  labs(x = "Number of curved setae/spines per pereiopod",
       y = "Standardised number of mates")

select_grads <- analyse_directional_selection %>% 
  filter(Sex == "Male") %>% 
  mutate(Selection_differential = z1_bar - z_bar) %>% 
  rename(Selection_gradient = beta_1) %>% 
  select(Species, Sex, Trait, contains("Selection")) %>% 
  pivot_longer(cols = contains("Selection"), values_to = "Strength", names_to = "Stat") %>% 
  mutate(Stat = factor(Stat, levels = c("Selection_gradient", "Selection_differential")),
         Trait = case_when(
           grepl("P1", Trait) ~ "Pereiopod 1",
           grepl("P2", Trait) ~ "Pereiopod 2",
           grepl("P3", Trait) ~ "Pereiopod 3",
           grepl("P4", Trait) ~ "Pereiopod 4",
           grepl("P5", Trait) ~ "Pereiopod 5",
           grepl("P6", Trait) ~ "Pereiopod 6",
           grepl("P7", Trait) ~ "Pereiopod 7",
           Trait == "Size" ~ "Size (in mm)",
           TRUE ~ Trait
         )) %>% 
  drop_na() %>% 
  ggplot(aes(x = Trait, y = Strength, fill = Species)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(values = c("forsmani" = "#3A9AB2", "praehirsuta" = "navy")) +
  facet_wrap(vars(Stat), scales = "free_y", ncol = 1) +
  labs(tag = "(B)") +
  theme_bw() +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank())



traits <- reproductive_success_parents %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species, P7.ep, P3.curv.setae),
            by = join_by("Parent_ID" == "ID_DNA_RAD", "Species")) %>% 
  pivot_longer(cols = c(Size, P7.ep, P3.curv.setae),
               values_to = "count", names_to = "Trait") %>% 
  filter(Sex == "Male") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, y = w1, colour = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.2) +
  geom_vline(data = analyse_directional_selection %>%
               filter(Sex == "Male",
                      Trait %in% c("Size", "P3.curv.setae", "P7.ep")) %>%
               select(Species, Sex, Trait, z1_bar, z_bar) %>%
               pivot_longer(contains("z"), names_to = "z_stat", values_to = "z_values"),
             aes(xintercept = z_values, colour = Species), lwd = 1.2) +
  geom_segment(data = analyse_directional_selection %>%
                 filter(Sex == "Male",
                        (Trait == "P7.ep") | (Trait == "Size" & Species == "forsmani")),
               aes(x = z_bar, xend = z1_bar, y = 2.25, colour = Species),
               arrow = arrow(length = unit(0.5, "cm")), lwd = 1.2) +
  scale_colour_manual(values = c("forsmani" = "#3A9AB2", "praehirsuta" = "navy")) +
  facet_grid2(Species ~ Trait, scales = "free",
              labeller = as_labeller(c(
                "forsmani" = "J. forsmani", "praehirsuta" = "J. praehirsuta",
                "P3.curv.setae" = "Pereiopod 3",
                "P7.ep" = "Pereiopod 7",
                "Size" = "Size (in mm)"
              ))) +
  labs(x = "Measure",
       y = "Relative mating success",
       tag = "(A)") +
  theme_bw() +
  theme(text = element_text(size = 30),
        strip.text.y = element_text(face = "italic", angle = 90),
        legend.position = "none")

ggarrange(traits, select_grads,
          nrow = 1,
          widths = c(1.5, 1))





reproductive_success_parents %>% 
  filter(Sex == "Male") %>% 
  drop_na() %>% 
  ggplot(aes(x = Size, y = w1)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female",
                      Trait == "Size"),
             aes(xintercept = z1_bar)) +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female",
                      Trait == "Size"),
             aes(xintercept = z_bar)) +
  geom_segment(data = analyse_directional_selection %>% 
                 filter(Sex != "Female",
                        Trait == "Size",
                        Species != "praehirsuta"),
               aes(x = z_bar, xend = z1_bar, y = 1.5),
               arrow = arrow(length = unit(0.3, "cm"))) +
  facet_grid2(Species ~ Trait, scales = "free") +
  my_theme +
  labs(x = "Size of individuals (mm)",
       y = "Relative mating success")


reproductive_success_parents %>% 
  filter(Sex == "Male") %>% 
  drop_na() %>% 
  ggplot(aes(x = Size, y = w1, colour = Species)) +
  geom_smooth(method = "lm", se = FALSE,
              lwd = 1.2) +
  geom_point(size = 4) +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female",
                      Trait == "Size"),
             aes(xintercept = z1_bar, colour = Species),
             lwd = 1.2) +
  geom_vline(data = analyse_directional_selection %>% 
               filter(Sex != "Female",
                      Trait == "Size"),
             aes(xintercept = z_bar, colour = Species),
             lwd = 1.2) +
  geom_segment(data = analyse_directional_selection %>% 
                 filter(Sex != "Female",
                        Trait == "Size",
                        Species == "forsmani"),
               aes(x = z_bar, xend = z1_bar, y = 2.6, colour = Species),
               arrow = arrow(length = unit(0.5, "cm")),
               lwd = 1.2) +
  scale_colour_manual(values = c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")) +
  labs(x = "Male phenotypic trait (e.g. number of setae)",
       y = "Relative mating success") +
  theme_bw() +
  theme(text = element_text(size = 30),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")


fm <- (reproductive_success_parents %>% 
         filter(Species == "forsmani",
                Sex == "Male") %>% 
         ggplot(aes(x = w1, y = w12)) +
         geom_point(size = 2, colour = "#3A9AB2") +
         geom_smooth(data = reproductive_success_parents %>% 
                       filter(Species == "forsmani", Sex == "Male", w1 > 0),
                     colour = "#3A9AB2", method = "lm", se = FALSE) +
         labs(tag = "(B)") +
         ylim(0, 4) +
         theme_bw() +
         theme(text = element_text(size = 20),
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.title = element_blank())) %>% 
  ggMarginal(type = "histogram", size = 10, fill = "#3A9AB2", bins = 20)

ff <- (reproductive_success_parents %>% 
         filter(Species == "forsmani",
                Sex == "Female") %>% 
         ggplot(aes(x = w1, y = w12)) +
         geom_point(size = 2, colour = "#3A9AB2") +
         geom_smooth(data = reproductive_success_parents %>% 
                       filter(Species == "forsmani", Sex == "Female", w1 > 0),
                     colour = "#3A9AB2", method = "lm", se = FALSE) +
         labs(y = "Standardised reproduction success",
              tag = "(A)") +
         xlim(0, 3) +
         ylim(0, 4) +
         theme_bw() +
         theme(text = element_text(size = 20),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank())) %>% 
  ggMarginal(type = "histogram", size = 10, fill = "#3A9AB2", bins = 20)

pm <- (reproductive_success_parents %>% 
         filter(Species == "praehirsuta",
                Sex == "Male") %>% 
         ggplot(aes(x = w1, y = w12)) +
         geom_point(size = 2, colour = "navy") +
         geom_smooth(data = reproductive_success_parents %>% 
                       filter(Species == "praehirsuta", Sex == "Male", w1 > 0),
                     colour = "navy", method = "lm", se = FALSE) +
         labs(x = "Standardised mating success",
              tag = "(D)") +
         ylim(0, 4) +
         theme_bw() +
         theme(text = element_text(size = 20),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())) %>% 
  ggMarginal(type = "histogram", size = 10, fill = "navy", bins = 20)

pf <- (reproductive_success_parents %>% 
         filter(Species == "praehirsuta",
                Sex == "Female") %>% 
         ggplot(aes(x = w1, y = w12)) +
         geom_point(size = 2, colour = "navy") +
         geom_smooth(data = reproductive_success_parents %>% 
                       filter(Species == "praehirsuta", Sex == "Female", w1 > 0),
                     colour = "navy", method = "lm", se = FALSE) +
         labs(x = "Standardised mating success",
              y = "Standardised reproduction success",
              tag = "(C)") +
         xlim(0, 3) +
         ylim(0, 4) +
         theme_bw() +
         theme(text = element_text(size = 20))) %>% 
  ggMarginal(type = "histogram", size = 10, fill = "navy", bins = 20)

ann1 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Females)"),
            parse = TRUE, size = 6, hjust = 0.1) +
  theme_void()
ann2 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Males)"),
            parse = TRUE, size = 6, hjust = 0.5) +
  theme_void()


ggarrange(ann1, ann2,
          ff, fm,
          pf, pm,
          ncol = 2, nrow = 3,
          widths = c(1, 0.88),
          heights = c(0.1, 0.88, 1)) %>% 
  ggsave(plot = ., "/shared/home/bpajot/Distribution_toto.png", scale =5, width = 1200, height = 1000, units = "px")


#### BIN ####
shapes <- shapes_cast() %>% 
  filter(set == "flower") %>% 
  pull(shape)

reproductive_success_parents %>% 
  ggplot(aes(x = w1, y = w12)) +
  geom_casting(aes(shape = factor(after_stat(n)), group = after_stat(n)),
               stat = "sum", size = 0.25, fill = "black") +
  geom_smooth(method = "lm", se = FALSE) +
  scale_shape_manual(values = shapes) +
  facet_grid2(Species ~ Sex, scales = "free") +
  my_theme +
  labs(shape = "Count")
