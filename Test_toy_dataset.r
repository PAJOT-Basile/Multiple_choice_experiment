########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "parallel", "ggh4x")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)


########################################################
################### Useful functions ###################
########################################################
"%not_in%" <- function(x, y){
  #' Not in
  #' 
  #' This function is the opposite of the "%in%" function. It checks if the 
  #' values contained in input x are NOT contained in y.
  #' 
  #'@param x (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the y vector
  #'@param y (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the x vector
  #'    
  #'@returns (bool)
  #'    This function returns TRUE for each value of x that is NOT in y and FALSE
  #'    for every value of x that is in y.
  #' @export
  
  return(!(x %in% y))
}

str_split_last <- function(string, pattern){
  #' Split string on last pattern
  #' 
  #' @description
  #' This function works like the `str_split` functions, that allows to cut
  #' a string on the last occurrence of a pattern in the string
  #' 
  #' @param string String to cut
  #' @param pattern Pattern on which to cut string
  #' 
  #' @return A string that contains all the provided string except what comes after
  #' the last occurrence of the pattern in the string.
  nb_substrings <- str_count(string, pattern)
  nb_substrings <- ifelse(pattern != "", nb_substrings + 1, nb_substrings)
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

compute_beta_and_S_in_bootstraps <- function(i, df, cols_to_take, onsets = c(1, 2, 3, 12, 123), resample = TRUE){
  if (resample){
    lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  }else{
    lines_to_sample <- 1:nrow(df)
  }
  df_samp <- df[lines_to_sample, ]
  
  # Calculate the parameters
  df_samp$p0 <- 1 / nrow(df_samp)
  df_samp$w1 <- ifelse(is.na(df_samp$W1 / mean(df_samp$W1, na.rm = TRUE)), 0, df_samp$W1 / mean(df_samp$W1, na.rm = TRUE))
  df_samp$p1 <- ifelse(is.na(df_samp$p0 * df_samp$w1), 0, df_samp$p0 * df_samp$w1)
  df_samp$w2 <- ifelse(is.na(df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE)), 0, df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE))
  df_samp$p2 <- ifelse(is.na(df_samp$p1 * df_samp$w2), 0, df_samp$p1 * df_samp$w2)
  df_samp$w3 <- ifelse(is.na(df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE)), 0, df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE))
  df_samp$p3 <- ifelse(is.na(df_samp$p2 * df_samp$w3), 0, df_samp$p2 * df_samp$w3)
  
  
  S_vals <- vector(mode = "list", length = length(cols_to_take))
  single_onsets <- onsets[sapply(onsets, function(x){length(str_split(x, "", simplify = TRUE))}, simplify = "vector") == 1] %>% 
    c(0)
  for (trait in cols_to_take){
    z_bars <- vector(mode = "list", length = length(single_onsets))
    names(z_bars) <- paste0("z", single_onsets, "_bar")
    for (onset in single_onsets){
      df_samp[[paste(trait, onset, sep = "_")]] <- df_samp[[trait]] * df_samp[[paste0("p", onset)]]
      z_bars[[which(onset == single_onsets)]] <- sum(df_samp[[paste(trait, onset, sep = "_")]], na.rm = TRUE)
      
    }
    
    S_vals[[which(trait == cols_to_take)]] <- data.frame(
      "Trait" = trait,
      "S1" = z_bars$z1_bar - z_bars$z0_bar,
      "S2" = z_bars$z2_bar - z_bars$z1_bar,
      "S3" = z_bars$z3_bar - z_bars$z2_bar,
      "S12" = z_bars$z2_bar - z_bars$z0_bar,
      "S123" = z_bars$z3_bar - z_bars$z0_bar
    )
  }
  S_vals <- S_vals %>% 
    bind_rows()
  
  names_beta_vect <- expand.grid(cols_to_take, onsets)
  beta_vect <- rep(NA, nrow(names_beta_vect))
  nb_traits <- length(cols_to_take)
  counter <- 0
  for (onset in onsets){
    corresponding_S_val <- paste0("S", onset)
    index_onset <- which(onset == onsets)
    onset_nb_1 <- case_when(onset %in% c(1, 2, 3) ~ onset - 1,
                            onset %in% c(12, 123) ~ 0)
    P_a <- df_samp %>% 
      select(all_of(paste0(cols_to_take, "_", onset_nb_1))) %>% 
      select_if(~sum(!is.na(.)) > 0)
    
    if(nrow(P_a) == 0){
      beta_vect[counter + 1:nb_traits] <- NA
    }else{
      P <- P_a %>% 
        cov(use = "pairwise") 
      
      index_not_null <- which(colSums(P) != 0)
      traits_to_keep <- colnames(P_a[index_not_null, index_not_null]) %>% 
        str_remove_all(paste0("_", onset_nb_1))
      # Use Tikhonov regularisation (add 1e-8 to diagonal of matrix, allowing to
      # inverse the P matrix all the time, even when it is nearly singular)
      P_inv <- (P[index_not_null, index_not_null]) %>% 
        solve()
      
      beta_vect[counter + index_not_null] <- P_inv %*% S_vals[S_vals$Trait %in% traits_to_keep, corresponding_S_val]
    }
    counter <- counter + nb_traits
  }
  names(beta_vect) <- paste0(names_beta_vect$Var1, "_", names_beta_vect$Var2)
  
  beta_vect %>%
    t() %>%
    as.data.frame() %>%
    mutate(Stat = "Beta") %>%
    rbind(S_vals %>%
            pivot_longer(-Trait, names_to = "Onset", values_to = "S") %>%
            mutate(Onset = paste(Trait, str_remove_all(Onset, "S"), sep = "_")) %>%
            select(-Trait) %>%
            pivot_wider(names_from = Onset, values_from = S) %>%
            mutate(Stat = "S")) %>%
    return()
}

bootstap_S_and_grads <- function(df, cols_to_take, nb_boots){
  nb_cores_to_use <- min(detectCores() - 1, nb_boots)
  cl <- makeCluster(nb_cores_to_use)
  clusterEvalQ(cl, library(tidyverse))
  
  df_out <- parLapply(cl = cl,
                      X = 1:nb_boots,
                      fun = compute_beta_and_S_in_bootstraps,
                      df, cols_to_take)
  stopCluster(cl)
  df_out <- bind_rows(df_out)
  return(df_out)
}

########################################################
################### Useful variables ###################
########################################################
set.seed(123456789)
nboots <- 10000

########################################################
###################### Create data #####################
########################################################
nb_samples <- 25

# Define the distribution parameters of the phenotypic traits
trait_parameters <- list(
  means = c("trait1" = 0,
            "trait2" = 2,
            "trait3" = 169),
  sds = c("trait1" = 1,
          "trait2" = 1,
          "trait3" = 1)
)

# Define the distribution parameters of the fitness parameters
fitness_parameters <- list(
  means = c("Episode1" = 1,
            "Episode2" = 18,
            "Episode3" = 10),
  sds = c("Episode1" = 0.5,
          "Episode2" = 0.5,
          "Episode3" = 2)
)


fitness_components <- lapply(c("trait", "Episode"), function(trait_or_episode, nb_samples, trait_params, fitness_params){
  lapply(c(1, 2, 3), function(trait_episode_nb, trait_or_episode, nb_samples, trait_params, fitness_params){
    
    if (trait_or_episode == "trait"){
      list_to_sample_from <- trait_params
    }else{
      list_to_sample_from <- fitness_params
    }
    name_sample_from <- paste0(trait_or_episode, trait_episode_nb)
    if (name_sample_from %in% c("trait2", "Episode2", "Episode1", "trait3")){
      trait_distribution <- rpois(nb_samples, lambda = list_to_sample_from[["means"]][name_sample_from])
    }else if (name_sample_from == "Episode3"){
      trait_distribution <- rbeta(nb_samples,
                                  shape1 = list_to_sample_from[["means"]][name_sample_from],
                                  shape2 = list_to_sample_from[["sds"]][name_sample_from])
    }else{
          trait_distribution <- rnorm(nb_samples,
                                mean = list_to_sample_from[["means"]][name_sample_from],
                                sd = list_to_sample_from[["sds"]][name_sample_from])
    }
    data.frame(
      "V1" = paste(trait_or_episode, trait_episode_nb, sep = "_"),
      "V2" = trait_distribution
    ) %>% 
      return()
  },
  trait_or_episode, nb_samples, trait_params, fitness_params) %>% 
    bind_rows()

}, nb_samples, trait_parameters, fitness_parameters) %>% 
  bind_rows() %>% 
  group_by(V1) %>% 
  mutate(V3 = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "V1", values_from = "V2") %>% 
  rename(ID_indiv = V3,
         W1 = Episode_1,
         W2 = Episode_2,
         W3 = Episode_3,
         traitA = trait_1,
         traitB = trait_2,
         traitC = trait_3) %>% 
  # Fitness cannot be negative so we change the values of W1 if they are negative to be equal to 0
  mutate(W1 = case_when(W1 < 0 ~ 0,
                        TRUE ~ W1),
  # And should the W1 value be changed to 0, individuals could not have any value for W2 or W3.
         across(c(W2, W3), ~ case_when(W1 == 0 ~ NA,
                                       . < 0 ~ 0,
                                      TRUE ~ .)))


########################################################
######### Selection differentials and gradients ########
########################################################
traits_to_use <- fitness_components %>% 
  select(starts_with("trait")) %>% 
  names() %>% 
  matrix(nrow = 1, ncol = 3, byrow = TRUE)


beta_and_S_vals <- apply(traits_to_use, 1, function(cols_to_take, df){
      
      compute_beta_and_S_in_bootstraps(1, df, cols_to_take, resample = FALSE) %>%
        pivot_longer(-Stat, names_to = "Trait", values_to = "Stat_value")
      
    }, fitness_components) %>%
      bind_rows() %>% 
  rowwise() %>%
  mutate(Onset = str_split_last(Trait, "_"),
         Trait = str_remove(Trait, "_[^_]*$")) %>% 
  pivot_wider(names_from = Stat, values_from = Stat_value)

########################################################
################# Confidence intervals #################
########################################################
distrib_sel_diffs_and_grads <- apply(traits_to_use, 1, function(cols_to_take, df, nb_boots){
      bootstap_S_and_grads(df, cols_to_take, nb_boots) %>%
        pivot_longer(-Stat, names_to = "Trait", values_to = "Stat_value")
      
    }, fitness_components, nboots) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(Onset = str_split_last(Trait, "_"),
         Trait = str_remove(Trait, "_[^_]*$")) %>%
  group_by(Trait, Onset, Stat) %>%
  mutate(n = row_number()) %>%
  pivot_wider(names_from = c(Stat, Onset), values_from = Stat_value) %>%
  ungroup() %>%
  select(-n)

CI_stats <- distrib_sel_diffs_and_grads %>% 
  pivot_longer(-Trait, names_to = "Stat", values_to = "Stat_value") %>% 
  group_by(Trait, Stat) %>% 
  summarize(CI_2.5 = quantile(Stat_value, probs = 0.025, na.rm = TRUE),
            CI_97.5 = quantile(Stat_value, probs = 0.975, na.rm = TRUE),
            mean_bootstrap_val = mean(Stat_value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Onset = str_split_fixed(Stat, "_", 2)[, 2],
         Stat = str_split_fixed(Stat, "_", 2)[, 1])


########################################################
################ Graphic representation ################
########################################################
stat_name_to_keep <- "Beta"


calc_val <- beta_and_S_vals %>% 
  pivot_longer(c(Beta, S), names_to = "Stat_name", values_to = "Stat") %>% 
  mutate(Stat_name = paste(Stat_name, Onset, sep = "_")) %>% 
  select(-Onset) %>% 
  filter(grepl(stat_name_to_keep, Stat_name))

# Get the CI
calc_ci <- CI_stats %>% 
  pivot_longer(c("CI_2.5", "CI_97.5"), names_to = "CI", values_to = "CI_val") %>% 
  mutate(Stat_name = paste(Stat, Onset, sep = "_")) %>% 
  filter(grepl(stat_name_to_keep, Stat_name))


distrib_sel_diffs_and_grads %>% 
  pivot_longer(-Trait, names_to = "Stat_name", values_to = "Distrib") %>% 
  filter(grepl(stat_name_to_keep, Stat_name)) %>% 
  ggplot(aes(x = Distrib)) +
  geom_histogram(bins = 100) +
  geom_vline(data = calc_val,
             aes(xintercept = Stat), colour = "red") +
  geom_vline(data = calc_ci,
             aes(xintercept = CI_val), colour = "green") +
  facet_grid2(Trait ~ Stat_name, scales = "free", independent = "x")



# Verification of the trait and fitness distributions
# Traits
fitness_components %>%
  select(starts_with("trait")) %>%
  pivot_longer(everything(),
               names_to = "Trait",
               values_to = "z") %>%
  ggplot(aes(x = z)) +
  geom_histogram(bins = 100) +
  facet_wrap(vars(Trait), ncol = 1, scales = "free")
# 
# fitness_components %>%
#   select(starts_with("trait")) %>% 
#   apply(2, sd, na.rm = TRUE)
# 
# fitness_components %>%
#   select(starts_with("trait")) %>% 
#   apply(2, mean, na.rm = TRUE)
# 
# # Fitness
# fitness_components %>%
#   select(starts_with("W")) %>%
#   mutate(W12 = ifelse(is.na(W2), 0, W1 * W2),
#          W123 = ifelse(is.na(W2), 0, W1 * W2 * W3)) %>% 
#   pivot_longer(everything(),
#                names_to = "Fitness",
#                values_to = "W") %>%
#   ggplot(aes(x = W)) +
#   geom_histogram(bins = 100) + 
#   facet_wrap(vars(Fitness), nrow = 1, scales = "free_y")
# 
# fitness_components %>%
#   select(starts_with("W")) %>%
#   mutate(W12 = ifelse(is.na(W2), 0, W1 * W2),
#          W123 = ifelse(is.na(W2), 0, W1 * W2 * W3)) %>% 
#   apply(2, sd, na.rm = TRUE)
# 
# fitness_components %>%
#   select(starts_with("W")) %>%
#   mutate(W12 = ifelse(is.na(W2), 0, W1 * W2),
#          W123 = ifelse(is.na(W2), 0, W1 * W2 * W3)) %>% 
#   apply(2, mean, na.rm = TRUE)
# 
# 
# diff_bootstraps_calc_value <- beta_and_S_vals %>% 
#   pivot_longer(c(Beta, S), names_to = "Stat", values_to = "calc_val") %>% 
#   left_join(CI_stats,
#             by = c("Trait", "Onset", "Stat")) %>% 
#   mutate(difference = calc_val - mean_bootstrap_val)
