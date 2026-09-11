########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "ggh4x", "parallel", "scales", "ggpubr", "stringi")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
###### Useful variables ######
##############################
set.seed(123456789)
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

indivs_to_remove <- c("E4961","E4706")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20))

theme_grid <- my_theme +
  theme(plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.margin = margin(t = 15, r = 15, b = 15, l = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15),
        text = element_text(size = 12))

nboots <- 50000

##############################
###### Useful functions ######
##############################
"%not_in%" <- function(x, y){return(!(x %in% y))}

compute_beta_and_S_in_bootstraps <- function(i, sex, species, df, cols_to_take, onsets = c(1, 2, 3, 5, 6), resample = TRUE){
  if (resample){
    # Sample randomly columns
    lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  }else{
    lines_to_sample <- 1:nrow(df)
  }
  df_samp <- df[lines_to_sample, ]
  
  # Calculate the parameters
  df_samp$w1 <- ifelse(is.na(df_samp$W1 / mean(df_samp$W1, na.rm = TRUE)), 0, df_samp$W1 / mean(df_samp$W1, na.rm = TRUE))
  df_samp$p1 <- ifelse(is.na(df_samp$p0 * df_samp$w1), 0, df_samp$p0 * df_samp$w1)
  df_samp$w2 <- ifelse(is.na(df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE)), 0, df_samp$W2 / sum(df_samp$p1 * df_samp$W2, na.rm = TRUE))
  df_samp$p2 <- ifelse(is.na(df_samp$p1 * df_samp$w2), 0, df_samp$p1 * df_samp$w2)
  df_samp$w3 <- ifelse(is.na(df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE)), 0, df_samp$W3 / sum(df_samp$p2 * df_samp$W3, na.rm = TRUE))
  df_samp$p3 <- ifelse(is.na(df_samp$p2 * df_samp$w3), 0, df_samp$p2 * df_samp$w3)
  df_samp$w5 <- ifelse(is.na(df_samp$W5 / sum(df_samp$p3 * df_samp$W5, na.rm = TRUE)), 0, df_samp$W5 / sum(df_samp$p3 * df_samp$W5, na.rm = TRUE))
  
  
  S_vals <- vector(mode = "list", length = length(cols_to_take))
  for (trait in cols_to_take){
    df_samp[[paste0(trait, "_0")]] <- df_samp$p0 * df_samp[[trait]]
    df_samp[[paste0(trait, "_1")]] <- df_samp$p1 * df_samp[[trait]]
    df_samp[[paste0(trait, "_2")]] <- df_samp$p2 * df_samp[[trait]]
    df_samp[[paste0(trait, "_3")]] <- df_samp$p3 * df_samp[[trait]]
    df_samp[[paste0(trait, "_5")]] <- df_samp$p3 * df_samp$w5 * df_samp[[trait]]
    
    z0_bar <- sum(df_samp[[paste0(trait, "_0")]], na.rm = TRUE)
    z1_bar <- sum(df_samp[[paste0(trait, "_1")]], na.rm = TRUE)
    z2_bar <- sum(df_samp[[paste0(trait, "_2")]], na.rm = TRUE)
    z3_bar <- sum(df_samp[[paste0(trait, "_3")]], na.rm = TRUE)
    z5_bar <- sum(df_samp[[paste0(trait, "_5")]], na.rm = TRUE)
    
    S_vals[[which(trait == cols_to_take)]] <- data.frame(
      "Trait" = trait,
      "S1" = z1_bar - z0_bar,
      "S2" = z2_bar - z1_bar,
      "S3" = z3_bar - z2_bar,
      "S5" = z5_bar - z3_bar,
      "Stot" = z5_bar - z0_bar
    )
  }
  S_vals <- S_vals %>% 
    bind_rows()
  
  names_beta_vect <- expand.grid(cols_to_take, onsets)
  beta_vect <- rep(NA, nrow(names_beta_vect))
  counter <- 0
  for (onset in onsets){
    index_onset <- which(onset == onsets)
    onset_nb_1 <- ifelse(onset == 5, 3, ifelse(onset > 5, 0, onset - 1))
    P_a <- df_samp %>% 
      select(all_of(paste0(cols_to_take, "_", onset_nb_1))) %>% 
      select_if(~sum(!is.na(.)) > 0)
    
    if(nrow(P_a) == 0){
      beta_vect[counter + c(1, 2, 3)] <- NA
    }else{
      P <- P_a %>% 
        cov(use = "pairwise") 
      
      index_not_null <- which(colSums(P) != 0)
      # Use Tikhonov regularisation (add 1e-7 to diagonal of matrix, allowing to
      # inverse the P matrix all the time, even when it is nearly singular)
      P_inv <- (P[index_not_null, index_not_null] + diag(1e-7, nrow(P))[index_not_null, index_not_null]) %>% 
        solve()
      
      beta_vect[counter + index_not_null] <- P_inv %*% S_vals[, index_onset + 1][index_not_null]
    }
    counter <- counter + 3
  }
  names(beta_vect) <- paste0(names_beta_vect$Var1, "_", names_beta_vect$Var2)
  
  beta_vect %>%
    t() %>% 
    as.data.frame() %>%
    mutate(Stat = "Beta") %>%
    rbind(S_vals %>%
            pivot_longer(-Trait, names_to = "Onset", values_to = "S") %>%
            mutate(Onset = paste(Trait, ifelse(grepl("tot", Onset), 6, str_remove_all(Onset, "S")), sep = "_")) %>%
            select(-Trait) %>%
            pivot_wider(names_from = Onset, values_from = S) %>%
            mutate(Stat = "S")) %>%
    return()
}

bootstap_S_and_grads <- function(sex, species, df, cols_to_take, nb_boots){
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(tidyverse))

  df_out <- parLapply(cl = cl,
                      X = 1:nb_boots,
                      fun = compute_beta_and_S_in_bootstraps,
                      sex, species, df, cols_to_take)
  stopCluster(cl)
  df_out <- bind_rows(df_out)
  return(df_out)
}

compute_covar_matrices <- function(onset_nb, species, sex, phenos){
  onset_nb_1 <- ifelse(onset_nb == 5, 3, onset_nb - 1)
  P <- phenos %>% 
    select(ends_with(paste0("_", onset_nb_1)))  %>% 
    select_if(~sum(!is.na(.)) > 0)
  P_a <- P %>% 
    cov(use = "pairwise")
  index_not_null <- unname(which(colSums(P_a) != 0))
  P_i <- P_a[index_not_null, index_not_null] 

  if (length(P_i) == 1){
    P_i <- P_i %>% 
      as.data.frame() %>% 
      mutate(Trait1 = names(which(colSums(P) != 0)),
             Trait2 = names(which(colSums(P) != 0)),
             Onset = onset_nb_1) %>% 
      rename(var_covar = ".")
  }else if(length(P_i) == 0){
    P_i <- data.frame(
      "Trait1" = NA,
      "Trait2" = NA,
      "var_covar" = NA,
      "Onset" = onset_nb_1
    )
  }else{
    P_i <- P_i %>% 
      as.data.frame() %>% 
      rownames_to_column("Trait1") %>% 
      pivot_longer(!Trait1, names_to = "Trait2", values_to = "var_covar") %>% 
      mutate(Onset = onset_nb_1)
  }
  return(P_i)
}

str_split_last <- function(string, pattern){
  nb_substrings <- str_count(string, pattern) + 1
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

plot_linerange_traits <- function(df, traits, stat_name, ylab = ""){
  forsmani <- df %>%
    filter(Species == "J. forsmani",
           Trait %in% traits,
           Stat == stat_name,
           !is.na(Signif_S)) %>% 
    ggplot() +
    geom_hline(aes(yintercept = 0), colour = "black", lwd = 1.2) +
    geom_linerange(aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5)) +
    geom_linerange(data = . %>% 
                     filter(Onset == "Total"),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) +
    facet_grid2(Sex + Trait ~ ., scales = "free", independent = "y", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE)) +
    facetted_pos_scales(y = list(
      TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
    )) +
    labs(y = ylab) +
    theme_grid  +
    theme(legend.text = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(margin = margin(r=2)),
          axis.ticks.length.y = unit(1, "mm"))
  
  praehirsuta <- df %>% 
    filter(Species == "J. praehirsuta",
           Trait %in% traits,
           Stat == stat_name,
           !is.na(Signif_S)) %>% 
    ggplot() +
    geom_hline(aes(yintercept = 0), colour = "black", lwd = 1.2) +
    geom_linerange(aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5)) +
    geom_linerange(data = . %>% 
                     filter(Onset == "Total"),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) +
    facet_grid2(Sex + Trait ~ ., scales = "free", independent = "y", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE)) +
    facetted_pos_scales(y = list(
      TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
    )) +
    labs(x = "Selection episode", y = ylab) +
    theme_grid  +
    theme(legend.text = element_text(face = "italic"),
          strip.background.x = element_blank(),
          strip.text.x = element_blank())
  
  ann1 <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "bold(italic(J.~forsmani))"),
              parse = TRUE, size = 6, angle = 90) +
    theme_void()
  ann2 <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "bold(italic(J.~praehirsuta))"),
              parse = TRUE, size = 6, angle = 90) +
    theme_void()
  
  p <- ggarrange(ann1, forsmani,
                 ann2, praehirsuta,
                 nrow = 2, ncol = 2,
                 widths = c(0.05, 1),
                 common.legend = TRUE,
                 legend = "right") 
  
  return(p)
}

plot_heatmap_var_covar <- function(df, traits){
  df_filt <- df %>% 
    filter(Trait1 %in% traits, Trait2 %in% traits)
  
  max_var <- max(df_filt$var_covar, na.rm = TRUE)
  min_var <- min(df_filt$var_covar, na.rm = TRUE)
  
  
  forsmani <- df_filt %>% 
    filter(Species == "forsmani") %>% 
    ggplot(aes(x = Trait1, y = Trait2, fill = var_covar)) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(name = "Variance\nCovariance",
                         colours = c("navy", "dodgerblue", "skyblue", "white", "orange2", "red", "firebrick"),
                         limits = c(min_var, max_var)) +
    facet_grid2(Sex ~ Onset, scales = "free", independent = "x", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE
    )) +
    theme_grid  +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank())
  
  praehirsuta <- df_filt %>% 
    filter(Species == "praehirsuta") %>% 
    ggplot(aes(x = Trait1, y = Trait2, fill = var_covar)) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(name = "Variance/\nCovariance",
                         colours = c("navy", "dodgerblue", "skyblue", "white", "orange2", "red", "firebrick"),
                         limits = c(min_var, max_var)) +
    facet_grid2(Sex ~ Onset, scales = "free", independent = "x", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE
    )) +
    facetted_pos_scales(x = list(Sex == "Female" ~ scale_x_discrete(labels = NULL))) +
    theme_grid  +
    theme(axis.title = element_blank(),
          strip.background.x = element_blank(),
          strip.text.x = element_blank())
  
  ann1 <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "bold(italic(J.~forsmani))"),
              parse = TRUE, size = 6, angle = 90) +
    theme_void()
  ann2 <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "bold(italic(J.~praehirsuta))"),
              parse = TRUE, size = 6, angle = 90) +
    theme_void()
  
  
  p <- ggarrange(ann1, forsmani,
                 ann2, praehirsuta,
                 nrow = 2, ncol = 2,
                 widths = c(0.05, 1),
                 common.legend = TRUE,
                 legend = "right")
  
  return(p)
}

##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the dulicates
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "parents") %>% 
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
         Species = factor(Species, levels = c("forsmani", "praehirsuta")),
         # Remove the ">" that are in the data
         across(starts_with(paste0("P", 1:7)), ~ str_remove_all(., ">") %>% trimws() %>% as.numeric())) %>% 
  unique() %>% 
  mutate(Sum_curved_setaeP1P5 = sum(P1.curv.setae, P2.curv.setae, P3.curv.setae, P4.curv.setae, P5.curv.setae, na.rm = TRUE),
         Sum_spinesP4P7 = sum(P4.small.ep, P5.ep, P6.ep, P7.ep, na.rm = TRUE)) %>% 
  # Remove two individuals:
  # Female E4961 produced one offspring that died before it was genotyped
  # Male E4706 died after 2 days spent in the aquarium. It was degraded so we were not able to have any phenotype information on this individual
  # and it did not produce any offspring
  filter(ID_DNA_RAD %not_in% indivs_to_remove) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, Sum_curved_setaeP1P5, Sum_spinesP4P7), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()


small_metadata <- metadata %>%
  select(ID_DNA_RAD, Sex, Species, starts_with("P"), Size, starts_with("Sum_"), starts_with("St_")) %>% 
  select(-c(Phenotype, Photo_ID)) %>% 
  mutate(Sex = ifelse(Sex == "F", "Female", "Male"))


phenotype_names <- small_metadata %>%
  select(starts_with(paste0("P", 1:7), ignore.case = FALSE), Size, starts_with("Sum"), starts_with("St_")) %>%
  names()

# Reproductive success
load("Associated_data/Rdata/reproductive_success_parents.Rdata")
reproductive_success_parents <- reproductive_success_parents %>% 
  rename(Sex = Parent) %>% 
  mutate(Sex = str_remove_all(Sex, "s")) %>% 
  right_join(small_metadata,
             by = join_by("Parent_ID" == "ID_DNA_RAD", "Sex", "Species")) %>% 
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
    p5 = p3 * w5,
    across(c(w1, w12, w123, w1234, w1235, p1, w2, p2, w3, p3, w4, w5, p5), ~ ifelse(is.na(.), 0, .)),
    across(all_of(phenotype_names), ~ . * p0, .names = "{col}_0"),
    across(all_of(phenotype_names), ~ . * p1, .names = "{col}_1"),
    across(all_of(phenotype_names), ~ . * p2, .names = "{col}_2"),
    across(all_of(phenotype_names), ~ . * p3, .names = "{col}_3"),
    across(all_of(phenotype_names), ~ . * p5, .names = "{col}_5"),
  ) %>% 
  ungroup()


######################################################
####### Selection differentials and gradients ########
######################################################
phenotypes_and_fitness <- reproductive_success_parents %>%
  select(-c(W12, W123, W1234, W1235)) %>%
  pivot_longer(-c(Parent_ID, Sex, Species, p0, p1, p2, p3, p5, starts_with("w")),
               names_to = "Trait", values_to = "z") %>%
  rowwise() %>%
  mutate(Onset = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_split_last(Trait, "_"), "default"),
         Trait = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_remove(Trait, "_[^_]*$"), Trait))

signif_S <- phenotypes_and_fitness %>%
  pivot_wider(names_from = "Onset", values_from = "z", names_prefix = "z") %>% 
  group_by(Sex, Species, Trait) %>% 
  summarize(missing = sum(is.na(z0)),
            count = n()) %>% 
  filter(missing != count) %>% 
  select(-c(missing, count)) %>% 
  left_join(phenotypes_and_fitness %>% 
              pivot_wider(names_from = "Onset", values_from = "z", names_prefix = "z"),
            by = c("Sex", "Species", "Trait")) %>% 
  group_by(Sex, Species, Trait) %>% 
  summarize(pval_S1 = t.test(z1, z0)$p.value,
            pval_S2 = t.test(z2, z1)$p.value,
            pval_S3 = t.test(z3, z2)$p.value,
            pval_S5 = t.test(z5, z3)$p.value,
            pval_STotal = t.test(z5, z0)$p.value) %>% 
  pivot_longer(starts_with("pval"), names_to = "Onset", values_to = "Signif_S") %>% 
  mutate(Onset = str_remove_all(Onset, "pval_S"),
         Signif_S = ifelse(Signif_S <= 0.05, "*", "")) 


traits_to_use <- matrix(c("Size", "Sum_curved_setaeP1P5", "Sum_spinesP4P7",
                          "St_Size", "St_Sum_curved_setaeP1P5", "St_Sum_spinesP4P7"),
                        nrow = 2, ncol = 3, byrow = TRUE)


beta_and_S_vals <- lapply(c("praehirsuta", "forsmani"), function(species, df, traits_to_use){
  df_filt <- df %>% 
    filter(Species == species)
  
  lapply(c("Male", "Female"), function(sex, species, df, traits_to_use){
    df_filt <- df %>% 
      filter(Sex == sex)
    
    apply(traits_to_use, 1, function(cols_to_take, sex, species, df){

      compute_beta_and_S_in_bootstraps(1, sex, species, df, cols_to_take, resample = FALSE) %>% 
        pivot_longer(-Stat, names_to = "Trait", values_to = "Stat_value")
      
    }, sex, species, df_filt) %>% 
      bind_rows() %>% 
      mutate(Sex = sex, Species = species)
    
  }, species, df_filt, traits_to_use) %>% 
    bind_rows()
  
}, reproductive_success_parents, traits_to_use) %>% 
  bind_rows() %>% 
  rowwise() %>%
  mutate(Onset = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_split_last(Trait, "_"), "default"),
         Onset = ifelse(Onset == "6", "Total", Onset),
         Trait = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_remove(Trait, "_[^_]*$"), Trait)) %>% 
  pivot_wider(names_from = Stat, values_from = Stat_value) %>% 
  left_join(signif_S, by = c("Sex", "Species", "Trait", "Onset"))

########################################
######### Confidence intervals #########
########################################
# distrib_sel_diffs_and_grads <-   lapply(c("praehirsuta", "forsmani"), function(species, df, traits_to_use, nb_boots){
#   print(paste0("Species: ", species))
#   df_filt <- df %>%
#     filter(Species == species)
# 
#   lapply(c("Male", "Female"), function(sex, species, df, traits_to_use, nb_boots){
#     print(paste0("  - Sex: ", sex))
#     df_filt <- df %>%
#       filter(Sex == sex)
# 
#     apply(traits_to_use, 1, function(cols_to_take, sex, species, df, nb_boots){
#       print(paste0("    - Traits: ", paste(cols_to_take, collapse = ", ")))
# 
#       bootstap_S_and_grads(sex, species, df, cols_to_take, nb_boots) %>%
#         pivot_longer(-Stat, names_to = "Trait", values_to = "Stat_value")
# 
#     }, sex, species, df_filt, nb_boots) %>%
#       bind_rows() %>%
#       mutate(Sex = sex, Species = species)
# 
#   }, species, df_filt, traits_to_use, nb_boots) %>%
#     bind_rows()
# 
# }, reproductive_success_parents, traits_to_use, nboots) %>%
#   bind_rows() %>%
#   rowwise() %>%
#   mutate(Onset = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_split_last(Trait, "_"), "default"),
#          Onset = ifelse(Onset == "6", "Total", Onset),
#          Trait = ifelse(grepl("_[0-9]$", Trait, perl = TRUE), str_remove(Trait, "_[^_]*$"), Trait)) %>%
#   group_by(Trait, Sex, Species, Onset, Stat) %>%
#   mutate(n = row_number()) %>%
#   pivot_wider(names_from = c(Stat, Onset), values_from = Stat_value) %>%
#   ungroup() %>%
#   select(-n)
# 
# 
# save(distrib_sel_diffs_and_grads, file = "Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000.rdata")
load("Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000.rdata")

CI_stats <- distrib_sel_diffs_and_grads %>% 
  pivot_longer(-c(Species, Sex, Trait), names_to = "Stat", values_to = "Stat_value") %>% 
  group_by(Species, Sex, Trait, Stat) %>% 
  summarize(CI_2.5 = quantile(Stat_value, probs = 0.025, na.rm = TRUE),
            CI_97.5 = quantile(Stat_value, probs = 0.975, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Onset = str_split_fixed(Stat, "_", 2)[, 2],
         Stat = str_split_fixed(Stat, "_", 2)[, 1])
########################################
############# Plot results #############
########################################
beta_and_S_vals <- beta_and_S_vals %>% 
  pivot_longer(c(Beta, S), names_to = "Stat", values_to = "Stat_value") %>% 
  left_join(CI_stats, by = c("Sex", "Species", "Trait", "Onset", "Stat")) %>% 
  mutate(Onset = case_when(Onset == "1" ~ "Relative number of\nreproductive partners (1)",
                           Onset == "2" ~ "Relative fecundity\n(2)",
                           Onset == "3" ~ "Relative contribution to brood\n(3)",
                           Onset == "5" ~ "Relative survival of offspring\nat 60 days (5)",
                           Onset == "Total" ~ "Total") %>% 
           factor(levels = c("Relative number of\nreproductive partners (1)", "Relative fecundity\n(2)",
                             "Relative contribution to brood\n(3)", "Relative survival of offspring\nat 60 days (5)", "Total")),
         Trait = case_when(Trait == "Size" ~ "Size\n(mm)",
                           Trait == "Sum_curved_setaeP1P5" ~ "Sum of curved\nsetae",
                           Trait == "Sum_spinesP4P7" ~ "Sum of spines",
                           Trait == "St_Size" ~ "Standardised size\n(mm)",
                           Trait == "St_Sum_curved_setaeP1P5" ~ "Standardised sum of curved\nsetae",
                           Trait == "St_Sum_spinesP4P7" ~ "Standardised sum of spines",) %>% 
           factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines", "Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines")),
         Species = paste0("J. ", Species) %>% factor(levels = c("J. forsmani", "J. praehirsuta")))


# Look at the selection differentials
selection_differentials_plot <- plot_linerange_traits(df = beta_and_S_vals,
                                                      traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines"),
                                                      stat_name = "S",
                                                      ylab = "Selection differentials")
selection_differentials_plot
ggsave(plot = selection_differentials_plot, "../Figures/tmp/Sexual_selection/Selection_differentials.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")


selection_differentials_plot_std <- plot_linerange_traits(df = beta_and_S_vals,
                                                          traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                          stat_name = "S",
                                                          ylab = "Selection differentials")
selection_differentials_plot_std
ggsave(plot = selection_differentials_plot_std, "../Figures/tmp/Sexual_selection/Selection_differentials_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

# Look at the selection gradients
selection_gradients_plot <- plot_linerange_traits(df = beta_and_S_vals,
                                                  traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines"),
                                                  stat_name = "Beta",
                                                  ylab = "Selection gradients")
selection_gradients_plot
ggsave(plot = selection_gradients_plot, "../Figures/tmp/Sexual_selection/Selection_gradients.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")


selection_gradients_plot_std <- plot_linerange_traits(df = beta_and_S_vals,
                                                          traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                          stat_name = "S",
                                                          ylab = "Selection differentials")
selection_gradients_plot_std
ggsave(plot = selection_gradients_plot_std, "../Figures/tmp/Sexual_selection/Selection_gradients_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")


########################################################
################# Covariance matrices ##################
########################################################
variance_covariance_matrices <- lapply(c("forsmani", "praehirsuta"), function(species, df){
  df_filt <- df %>% 
    filter(Species == species)
  lapply(c("Male", "Female"), function(sex, species, df){
    df_filt <- df %>% 
      filter(Sex == sex)
    lapply(c(1, 2, 3, 5), compute_covar_matrices, species, sex, df_filt) %>% 
      bind_rows() %>% 
      mutate(Sex = sex,
             Species = species)
  }, species, df_filt) %>% 
    bind_rows()
}, reproductive_success_parents %>% 
  select(Parent_ID, Sex, Species, starts_with("Sum"), starts_with("St_S"), starts_with("Size"))) %>% 
  bind_rows() %>% 
  mutate(across(starts_with("Trait"), ~ case_when(stri_replace_last_regex(., "_[0-9]", "") == "Size" ~ "Size\n(mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_curved_setaeP1P5" ~ "Sum of curved setae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_spinesP4P7" ~ "Sum of spines",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Size" ~ "Standardised\nsize (mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_curved_setaeP1P5" ~ "Standardised sum\nof curved setae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_spinesP4P7" ~ "Standardised sum\nof spines",
                                                  TRUE ~ NA)),
         Onset = case_when(Onset == 0 ~ "Relative number of\nreproductive partners (1)",
                           Onset == 1 ~ "Relative fecundity (2)",
                           Onset == 2 ~ "Relative contribution\nto brood\n(3)",
                           Onset == 3 ~ "Relative survival\nof offspring at 60 days (5)") %>% 
           factor(levels = c("Relative number of\nreproductive partners (1)", "Relative fecundity (2)", "Relative contribution\nto brood\n(3)", "Relative survival\nof offspring at 60 days (5)")))



var_covar_plot <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                       traits = c("Size\n(mm)", "Sum of curved setae", "Sum of spines"))
var_covar_plot
ggsave(plot = var_covar_plot, "../Figures/tmp/Sexual_selection/Variance_covariance.png",
         grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

var_covar_plot_std <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                                             traits = c("Standardised\nsize (mm)", "Standardised sum\nof curved setae", "Standardised sum\nof spines"))
var_covar_plot_std
ggsave(plot = var_covar_plot_std, "../Figures/tmp/Sexual_selection/Variance_covariance_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

#####################################################
########## Plots specific cases gradients ###########
#####################################################

plot_phenotype_fitness_regression <- function(df, phenotype, fitness){
  fitness_onset <- c("b4" = ifelse(fitness == "tot", 1, ifelse(fitness == "w5", 3, as.numeric(str_split_fixed(fitness, "", 2)[, 2]) - 1)),
                     "after" = ifelse(fitness == "tot", 5, as.numeric(str_split_fixed(fitness, "", 2)[, 2])))
  fitness <- ifelse(fitness == "tot", "w5", fitness)
  
  phenotypes <- c("b4" = paste0(phenotype, "_", fitness_onset["b4"]), "after" = paste0(phenotype, "_", fitness_onset["after"]))
  
  start_arrow <- sum(df[[phenotypes["b4"]]], na.rm = TRUE)
  end_arrow <- sum(df[[phenotypes["after"]]], na.rm = TRUE)
  height_arrow <- 1.5*mean(df[[fitness]], na.rm = TRUE)
  
  sum_lm <- lm(as.formula(paste(fitness, "~", phenotype)), data = df) %>% 
    summary()
  
  
  r_squared <- sum_lm$r.squared
  beta_reg <- sum_lm$coefficients[2, 1]
  intercept <- sum_lm$coefficients[1, 1]
  p_val <- sum_lm$coefficients[2, 4]

  if(intercept < 0){
    fmt <- "atop(R^2 == %.2f *','~ p == %.2e , %s == %.2f*z~%.2f)"
  }else{
    fmt <- "atop(R^2 == %.2f *','~ p == %.2e , %s == %.2f*z~+~%.2f)"
  }
  lab <- sprintf(fmt, r_squared, p_val, fitness, beta_reg, intercept)
  
  
  p <- df %>% 
    ggplot() +
    geom_point(aes(x = !!sym(phenotype), y = !!sym(fitness)), colour = "blue", size = 3.5, alpha = 0.5) +
    geom_vline(aes(xintercept = sum(!!sym(phenotypes["b4"]), na.rm = TRUE)), colour = "blue") +
    geom_vline(aes(xintercept = sum(!!sym(phenotypes["after"]), na.rm = TRUE)), colour = "red") +
    annotate("segment", x = start_arrow, xend = end_arrow, y = height_arrow,
             arrow = arrow(length = unit(0.5, "cm"))) +
    annotate("text", x = Inf, y = Inf, label = lab,
             hjust = 1.2, vjust = 1.5, size = 7, parse = TRUE) +
    geom_smooth(aes(x = !!sym(phenotype), y = !!sym(fitness)), method = "lm", se = FALSE, colour = "black") +
    my_theme
  return(p)
}


save_plot_regression <- function(df, sex, species, pheno, fitness){
  x_text <- case_when(pheno == "Size" ~ "Size (mm)",
                      pheno == "Sum_spinesP4P7" ~ "Total number of spines",
                      pheno == "Sum_curved_setaeP1P5" ~ "Total number of curved setae")
  y_text <- case_when(fitness == "w1" ~ "Relative number of\nreproductive partners (1)",
                      fitness == "w2" ~ "Relative fecundity (2)",
                      fitness == "w3" ~ "Relative contribution to brood\n(3)",
                      fitness == "w5" ~ "Relative survival of offspring\nat 60 days (5)",
                      fitness == "tot" ~ "Total relative fitness\n(Relative number of offspring alive at 60 days)")
  
  p <- df %>% 
    filter(Sex == sex, Species == species) %>% 
    plot_phenotype_fitness_regression(pheno, fitness) +
    labs(
      x = x_text,
      caption = paste(sex, species))
  ggsave(filename = paste0("../Figures/Suppfigs/tmp/", paste(pheno, fitness, sex,  species, sep = "_"), ".png"),
         plot = p,
         device = grDevices::png,
         width = 1000, height = 700,
         units = "px", scale = 4)
  
  return(p)
}


save_plot_regression(reproductive_success_parents, "Male", "praehirsuta", "Sum_spinesP4P7", "tot")

reproductive_success_parents %>% 
  filter(Sex == "Female", Species == "praehirsuta") %>% 
  plot_phenotype_fitness_regression("Size", "tot")

######
reproductive_success_parents %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  ggplot(aes(x = Sum_curved_setaeP1P5, y = Size)) +
  geom_point(aes(colour = w1), size = 5, alpha = 0.7) +
  scale_colour_gradientn(colours = c("#A4A4A4", "orange2", "red", "firebrick")) +
  my_theme

toto <- reproductive_success_parents %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  mutate(toto = case_when(w1235 < 0.5 ~ "#A4A4A4",
                          w1235 >= 0.5 & w1235 < 1 ~ "orange2",
                          w1235 >= 1 & w1235 < 1.5 ~ "red",
                          w1235 > 1.5 ~ "firebrick"))

plot3d(x = toto$Sum_curved_setaeP1P5, y = toto$Size, z = toto$w1235,
       type = "s", col = toto$toto)
rglwidget(width = 520, height = 520)
