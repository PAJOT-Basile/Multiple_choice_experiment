# Decompose the opportunity for selection over three episodes of selection

########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "parallel")
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

nboots <- 50000

########################################
########### Useful functions ###########
########################################
"%not_in%" <- function(x, y) return(!(x %in% y))

compute_variances_in_bootstrap <- function(x, df){
  # Sample randomly columns
  lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  df_samp <- df[lines_to_sample, ]
  
  # Calculate the parameters
  df_samp$w1 <- df_samp$W1 / mean(df_samp$W1, na.rm = TRUE)
  df_samp$w12 <- df_samp$W12 / mean(df_samp$W12, na.rm = TRUE)
  df_samp$w123 <- df_samp$W123 / mean(df_samp$W123, na.rm = TRUE)
  
  df_samp$p1 <- df_samp$p0 * df_samp$w1
  df_samp$w2 <- df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE)
  
  df_samp$p2 <- df_samp$p1 * df_samp$w2
  df_samp$w3 <- df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE)
  
  df_samp$p3 <- df_samp$p2 * df_samp$w3
  df_samp <- df_samp %>% 
    mutate(across(c(w1, w12, w123, p1, w2, p2, w3, p3), ~ ifelse(is.na(.), 0, .)))

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
  
  to_ret <- c("I1" = I1, "I2" = I2, "I12" = I12,
              "COI1.2" = COI1.2, "COI1.2_1" = COI1.2_1, "COI12.2_1" = COI12.2_1, "COI12.2" = COI12.2,
              "I3" = I3, "I123" = I123,
              "COI1.3" = COI1.3, "COI2.3" = COI2.3, "COI12.3" = COI12.3, "COI12.3_2" = COI12.3_2, "COI123.3_2" = COI123.3_2, "COI123.3" = COI123.3)
  
  to_ret %>%
    t() %>% 
    as_tibble() %>% 
  return()
}

bootstrap_variances <- function(df, nb_boots){
  cl <- makeCluster(detectCores() - 1)
  clusterEvalQ(cl, library(tidyverse))
  df_out <- parLapply(cl = cl,
            X = 1:nb_boots,
            fun = compute_variances_in_bootstrap,
            df) %>% 
    bind_rows() 
  
  stopCluster(cl)
  return(df_out)
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
         across(c(contains("I", ignore.case = FALSE), contains("C")), ~ round(., digits = 2)))

########################################
######### Confidence intervals #########
########################################
# distrib_opportunity_selection <- lapply(c("forsmani", "praehirsuta"), function(species){
#   lapply(c("M", "F"), function(sex, sp){
#     fitness_components %>%
#       filter(Sex == sex, Species == sp) %>%
#       bootstrap_variances(., nboots) %>%
#       mutate(Species = sp, Sex = sex)
#   }, species)
# }) %>%
#   bind_rows() %>%
#   mutate(C1 = COI1.2 + COI1.2_1 + (COI12.2_1 - COI12.2),
#          C2 = COI12.3 + COI12.3_2 + (COI123.3_2 - COI123.3))
# 
# save(distrib_opportunity_selection, file = "Associated_data/Rdata/distribution_variances_50000_boots_three_episodes.rda")
load("Associated_data/Rdata/distribution_variances_50000_boots_three_episodes.rda")


CI_stats <- distrib_opportunity_selection %>% 
  pivot_longer(cols = !c(Sex, Species), names_to = "Stat_name", values_to = "Stat") %>% 
  group_by(Species, Sex, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat, probs = 0.025),
            CI_97.5 = quantile(Stat, probs = 0.975))


#########################################################
######### Look at the values to fill up table 2 #########
#########################################################
I_vals %>% 
  select(-starts_with("Percentage")) %>% 
  pivot_longer(!c(Species, Sex), names_to = "Stat_name", values_to = "Stat") %>% 
  left_join(CI_stats,
            by = c("Species", "Sex", "Stat_name")) %>% 
  left_join(I_vals %>%
              select(Species, Sex, starts_with("Percentage")) %>%
              pivot_longer(!c(Species, Sex),
                           names_to = "Stat_name",
                           values_to = "Percentage") %>%
              mutate(Stat_name = str_split_fixed(Stat_name, "_", 2)[, 2]),
            by = c("Species", "Sex", "Stat_name")) %>% 
  rowwise() %>% 
  mutate(CI = paste0("[", round(CI_2.5, digits = 2), "; ", round(CI_97.5, digits = 2), "]")) %>% 
  select(-c(CI_2.5, CI_97.5)) %>% 
  relocate(Percentage, .after = CI) %>% 
  filter(Stat_name == "I2") 



