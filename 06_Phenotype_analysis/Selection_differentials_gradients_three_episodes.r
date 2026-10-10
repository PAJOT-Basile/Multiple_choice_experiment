########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "ggh4x", "parallel", "scales", "ggpubr", "stringi")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################################
################### Useful variables ###################
########################################################
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

########################################################
################### Useful functions ###################
########################################################
source("98_Transversal_scripts/Transversal_functions.r")

# Weighted covariance of raw traits with weights p (the p_k of the onset)
weighted_cov <- function(Z, p){
  ok <- p > 0 & complete.cases(Z)
  if (sum(ok) < 2) return(NULL)
  Z <- Z[ok, , drop = FALSE]
  p <- p[ok] / sum(p[ok])
  mu <- colSums(Z * p)
  Zc <- sweep(Z, 2, mu)
  crossprod(Zc * sqrt(p))          # t(Zc) %*% diag(p) %*% Zc
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
    onset_nb_1 <- case_when(onset %in% c(1, 2, 3) ~ onset - 1,
                            onset %in% c(12, 123) ~ 0)
    
    beta <- rep(NA, nb_traits)
    # Get the columns of interest
    Z <- df_samp %>% 
      select(all_of(cols_to_take)) %>%
      select_if(~ var(., na.rm = TRUE) != 0) %>% 
      as.matrix()
    
    p <- df_samp %>% 
      select(paste0("p", onset_nb_1)) %>% 
      as.matrix()
    P <- weighted_cov(Z, p)
      
    keep <- which(diag(P) > 0)
    
    if(nrow(Z) == 0 | length(keep) == 0){
      beta_vect[counter + 1:nb_traits] <- beta
    }else{
      Pk <- P[keep, keep, drop = FALSE]
      
      ev <- eigen(Pk, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) < 1e-8 * max(ev)){
        beta_vect[counter + 1:nb_traits] <- beta   # relative, scale-free test
      } else {
        traits_to_keep <- colnames(Pk)
        
        res <- tryCatch(as.vector(solve(Pk, S_vals[S_vals$Trait %in% traits_to_keep, corresponding_S_val])),
                        error = function(e) NULL)
        
        if (is.null(res) || any(!is.finite(res))){
          beta_vect[counter + 1:nb_traits] <- NA
        } else {
          beta_vect[counter + keep] <- res
        }
      }
    counter <- counter + nb_traits
    }
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

bootstap_S_and_grads <- function(sex, species, df, cols_to_take, nb_boots){
  nb_cores_to_use <- min(detectCores() - 1, nb_boots)
  cl <- makeCluster(nb_cores_to_use)
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, varlist = c("weighted_cov"))
  
  df_out <- parLapply(cl = cl,
                      X = 1:nb_boots,
                      fun = compute_beta_and_S_in_bootstraps,
                      df, cols_to_take)
  stopCluster(cl)
  df_out <- bind_rows(df_out)
  df_out$Sex <- sex
  df_out$Species <- species
  return(df_out)
}

compute_covar_matrices <- function(onset_nb, species, sex, phenos){
  onset_nb_1 <- case_when(onset_nb %in% c(1, 2, 3) ~ onset_nb - 1,
                          onset_nb %in% c(12, 123) ~ 0)
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
                     filter(Onset %in% c("12\nFecundity", "123\nReproductive success")),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1.2) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) +
    facet_grid2(Sex + Trait ~ ., scales = "free", independent = "y", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE)) +
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
                     filter(Onset %in% c("12\nFecundity", "123\nReproductive success")),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1.2) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) +
    facet_grid2(Sex + Trait ~ ., scales = "free", independent = "y", strip = strip_nested(
      background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
      by_layer_y = TRUE)) +
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

plot_phenotype_fitness_regression <- function(df, phenotype, fitness){
  fitness_onset <- c("b4" = ifelse(fitness %in% c("w12", "w123"), 0, as.numeric(str_split_fixed(fitness, "", 2)[, 2]) - 1),
                     "after" = as.numeric(str_split_last(fitness, "")))

  phenotype <- ifelse(grepl("setae", phenotype), paste0(phenotype, "_default"), phenotype)
  
  phenotypes <- c("b4" = paste0(phenotype, "_", fitness_onset["b4"]),
                  "after" = paste0(phenotype, "_", fitness_onset["after"]))
  
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
  y_text <- case_when(fitness == "w1" ~ "1\nMating success",
                      fitness == "w2" ~ "2\nFecundity per mate",
                      fitness == "w12" ~ "12\nFecundity",
                      fitness == "w3" ~ "3\nOffspring survival",
                      fitness == "w123" ~ "123\nReproductive success") %>% 
    factor(levels = c("1\nMating success", "2\nFecundity per mate",
                      "12\nFecundity", "3\nOffspring survival", "123\nReproductive success"))
  x_text <- case_when(pheno == "Size" ~ "Size\n(mm)",
                      pheno == "Sum_curved_setaeP1P5" ~ "Sum of curved\nsetae",
                      pheno == "Sum_spinesP6P7" ~ "Sum of spines",
                      pheno == "St_Size" ~ "Standardised size\n(mm)",
                      pheno == "St_Sum_curved_setaeP1P5" ~ "Standardised sum of curved\nsetae",
                      pheno == "St_Sum_spinesP6P7" ~ "Standardised sum of spines") %>% 
    factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines", "Standardised size\n(mm)",
                      "Standardised sum of curved\nsetae", "Standardised sum of spines"))
  
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

########################################################
###################### Import data #####################
########################################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "parents") %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b")) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()


small_metadata <- metadata %>%
  select(ID_DNA_RAD, Sex, Species, starts_with("P"), starts_with("default"), Size, starts_with("Sum_"), starts_with("St_")) %>% 
  select(-c(Phenotype, Photo_ID, contains("offspring"))) %>%
  mutate(Sex = ifelse(Sex == "F", "Female", "Male"))


phenotype_names <- small_metadata %>%
  select(starts_with("trunc"), starts_with("P"), starts_with("default"), Size, starts_with("Sum"), starts_with("St_")) %>%
  names()

# Reproductive success
load("Associated_data/Rdata/Fitness_components_parents.rdata")
fitness_components <- fitness_components %>% 
  mutate(Sex = ifelse(Sex == "M", "Male", "Female")) %>% 
  right_join(small_metadata,
             by = join_by("Parent_ID" == "ID_DNA_RAD", "Sex", "Species")) %>% 
  group_by(Species, Sex) %>% 
  rename(W1 = nb_partners,
         W12 = nb_offspring,
         W123 = nb_offspring_surv60) %>% 
  select(-c(max_nb_offspring, nb_offspring_surv40, -contains("trunc"))) %>% 
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
    across(c(w1, w12, w123, p1, w2, p2, w3, p3), ~ ifelse(is.na(.), 0, .)),
    across(all_of(phenotype_names), ~ . * p0, .names = "{col}_0"),
    across(all_of(phenotype_names), ~ . * p1, .names = "{col}_1"),
    across(all_of(phenotype_names), ~ . * p2, .names = "{col}_2"),
    across(all_of(phenotype_names), ~ . * p3, .names = "{col}_3"),
  ) %>% 
  ungroup()


########################################################
######### Selection differentials and gradients ########
########################################################
phenotypes_and_fitness <- fitness_components %>%
  select(-c(W12, W123)) %>%
  pivot_longer(-c(Parent_ID, Sex, Species, p0, p1, p2, p3, starts_with("w")),
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
            pval_S12 = t.test(z2, z0)$p.value,
            pval_S123 = t.test(z3, z0)$p.value) %>% 
  pivot_longer(starts_with("pval"), names_to = "Onset", values_to = "Signif_S") %>% 
  mutate(Onset = str_remove_all(Onset, "pval_S"),
         pval = Signif_S,
         Signif_S = ifelse(Signif_S <= 0.05, "*", "")) 


traits_to_use <- matrix(c("Size", "Sum_curved_setaeP1P5_default", "Sum_spinesP6P7",
                          "St_Size", "St_Sum_curved_setaeP1P5_default", "St_Sum_spinesP6P7"),
                        nrow = 2, ncol = 3, byrow = TRUE)


beta_and_S_vals <- lapply(c("praehirsuta", "forsmani"), function(species, df, traits_to_use){
  df_filt <- df %>% 
    filter(Species == species)
  
  lapply(c("Male", "Female"), function(sex, species, df, traits_to_use){
    df_filt <- df %>% 
      filter(Sex == sex)

    apply(traits_to_use, 1, function(cols_to_take, sex, species, df){

      compute_beta_and_S_in_bootstraps(1, df, cols_to_take, resample = FALSE) %>%
        pivot_longer(-Stat, names_to = "Trait", values_to = "Stat_value")

    }, sex, species, df_filt) %>%
      bind_rows() %>%
      mutate(Sex = sex, Species = species)
    
  }, species, df_filt, traits_to_use) %>% 
    bind_rows()
  
}, fitness_components %>% 
  select(-c(ends_with(paste0("_", 0:3)), starts_with("p", ignore.case = FALSE), starts_with("w", ignore.case = FALSE))),
traits_to_use) %>% 
  bind_rows() %>% 
  rowwise() %>%
  mutate(Onset = str_split_last(Trait, "_"),
         Trait = str_remove(Trait, "_[^_]*$")) %>% 
  pivot_wider(names_from = Stat, values_from = Stat_value) %>% 
  left_join(signif_S, by = c("Sex", "Species", "Trait", "Onset"))

########################################################
################# Confidence intervals #################
########################################################
# distrib_sel_diffs_and_grads <- lapply(c("praehirsuta", "forsmani"), function(species, df, traits_to_use, nb_boots){
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
#         pivot_longer(-c(Stat, Sex, Species), names_to = "Trait", values_to = "Stat_value")
# 
#     }, sex, species, df_filt, nb_boots) %>%
#       bind_rows()
# 
#   }, species, df_filt, traits_to_use, nb_boots) %>%
#     bind_rows()
# 
# }, fitness_components %>% 
#   select(-c(ends_with(paste0("_", 0:3)), starts_with("p", ignore.case = FALSE), starts_with("w", ignore.case = FALSE))),
# traits_to_use, nboots) %>%
#   bind_rows() %>%
#   rowwise() %>%
#   mutate(Onset = str_split_last(Trait, "_"),
#          Trait = str_remove(Trait, "_[^_]*$")) %>%
#   group_by(Trait, Sex, Species, Onset, Stat) %>%
#   mutate(n = row_number()) %>%
#   pivot_wider(names_from = c(Stat, Onset), values_from = Stat_value) %>%
#   ungroup() %>%
#   select(-n)
# 
# 
# save(distrib_sel_diffs_and_grads, file = "Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000_three_episodes.rdata")
load("Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000_three_episodes.rdata")

CI_stats <- distrib_sel_diffs_and_grads %>% 
  pivot_longer(-c(Species, Sex, Trait), names_to = "Stat", values_to = "Stat_value") %>% 
  group_by(Species, Sex, Trait, Stat) %>% 
  summarize(CI_2.5 = quantile(Stat_value, probs = 0.025, na.rm = TRUE),
            CI_97.5 = quantile(Stat_value, probs = 0.975, na.rm = TRUE),
            mean_val_boots = mean(Stat_value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Onset = str_split_fixed(Stat, "_", 2)[, 2],
         Stat = str_split_fixed(Stat, "_", 2)[, 1])


stat_name <- "Beta"

calc_val <- beta_and_S_vals %>% 
  filter(Sex == "Male", Species == "forsmani", !grepl("St_", Trait)) %>% 
  pivot_longer(c(Beta, S), names_to = "Stat", values_to = "Stat_value") %>% 
  filter(grepl(stat_name, Stat)) %>% 
  mutate(Species = str_remove_all(Species, "J. "),
         Trait = case_when(grepl("Size", Trait) ~ "Size",
                           grepl("setae", Trait) ~ "Sum_curved_setaeP1P5_default",
                           TRUE ~ "Sum_spinesP6P7"),
         Onset = as.numeric(gsub("([0-9]+).*$", "\\1", Onset)),
         Stat_name = paste(Stat, Onset, sep = "_"))# %>% 
  # filter(Stat_name == "Beta_2") 

calc_ci <- CI_stats %>%
  filter(Sex == "Male", Species == "forsmani", !grepl("St_", Trait)) %>% 
  pivot_longer(contains("CI_"), names_to = "CI", values_to = "CI_val") %>% 
  filter(grepl(stat_name, Stat)) %>% 
  mutate(Stat_name = paste(Stat, Onset, sep = "_"))

distrib_sel_diffs_and_grads %>% 
  filter(Sex == "Male", Species == "forsmani", !grepl("St_", Trait)) %>% 
  pivot_longer(-c(Sex, Species, Trait), names_to = "Stat_name", values_to = "Stat") %>% 
  filter(grepl(stat_name, Stat_name)) %>% 
  # filter(Stat_name == "Beta_2") %>% 
  ggplot(aes(x = Stat)) + 
  geom_histogram(bins = 100) +
  geom_vline(data = calc_val,
             aes(xintercept = Stat_value),
             colour = "red") +
  geom_vline(data = calc_ci,
             aes(xintercept = CI_val),
             colour = "green") +
  facet_grid2(Trait ~ Stat_name, scales = "free", independent = "all")

########################################################
##################### Plot results #####################
########################################################
beta_and_S_vals <- beta_and_S_vals %>%
  pivot_longer(c(Beta, S), names_to = "Stat", values_to = "Stat_value") %>% 
  left_join(CI_stats, by = c("Sex", "Species", "Trait", "Onset", "Stat")) %>% 
  mutate(Onset = case_when(Onset == "1" ~ "1\nMating success",
                           Onset == "2" ~ "2\nFecundity per mate",
                           Onset == "12" ~ "12\nFecundity",
                           Onset == "3" ~ "3\nOffspring survival",
                           Onset == "123" ~ "123\nReproductive success") %>% 
           factor(levels = c("1\nMating success", "2\nFecundity per mate",
                             "12\nFecundity", "3\nOffspring survival", "123\nReproductive success")),
         Trait = case_when(Trait == "Size" ~ "Size\n(mm)",
                           Trait == "Sum_curved_setaeP1P5_default" ~ "Sum of curved\nsetae",
                           Trait == "Sum_spinesP6P7" ~ "Sum of spines",
                           Trait == "St_Size" ~ "Standardised size\n(mm)",
                           Trait == "St_Sum_curved_setaeP1P5_default" ~ "Standardised sum of curved\nsetae",
                           Trait == "St_Sum_spinesP4P7" ~ "Standardised sum of spines") %>% 
           factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines", "Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines")),
         Species = paste0("J. ", Species) %>% factor(levels = c("J. forsmani", "J. praehirsuta")))


# Look at the selection differentials
selection_differentials_plot <- plot_linerange_traits(df = beta_and_S_vals,
                                                      traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines"),
                                                      stat_name = "S",
                                                      ylab = "Selection differentials")

selection_differentials_plot
ggsave(plot = selection_differentials_plot, "../Figures/Selection_differentials.png",
       grDevices::png, scale = 3, width = 1200, height = 1300, units = "px")


selection_gradients_plot <- plot_linerange_traits(df = beta_and_S_vals,
                                                  traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines"),
                                                  stat_name = "Beta",
                                                  ylab = "Selection gradients")
selection_gradients_plot
ggsave(plot = selection_gradients_plot, "../Figures/Selection_gradients.png",
       grDevices::png, scale = 3, width = 1200, height = 1300, units = "px")




# Look at the standardised traits
selection_differentials_plot_std <- plot_linerange_traits(df = beta_and_S_vals,
                                                          traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                          stat_name = "S",
                                                          ylab = "Selection differentials")
selection_differentials_plot_std
ggsave(plot = selection_differentials_plot_std, "../Figures/tmp/Sexual_selection/Selection_differentials_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

selection_gradients_plot_std <- plot_linerange_traits(df = beta_and_S_vals,
                                                      traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                      stat_name = "S",
                                                      ylab = "Selection differentials")
selection_gradients_plot_std
ggsave(plot = selection_gradients_plot_std, "../Figures/tmp/Sexual_selection/Selection_gradients_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")


########################################################
################## Covariance matrices ##################
########################################################
variance_covariance_matrices <- lapply(c("forsmani", "praehirsuta"), function(species, df){
  df_filt <- df %>% 
    filter(Species == species)
  lapply(c("Male", "Female"), function(sex, species, df){
    df_filt <- df %>% 
      filter(Sex == sex)
    lapply(c(1, 2, 3), compute_covar_matrices, species, sex, df_filt) %>% 
      bind_rows() %>% 
      mutate(Sex = sex,
             Species = species)
  }, species, df_filt) %>% 
    bind_rows()
}, fitness_components %>% 
  select(Parent_ID, Sex, Species, starts_with("Sum"), starts_with("St_S"), starts_with("Size")) %>% 
  select(-contains("trunc"))) %>% 
  bind_rows() %>% 
  mutate(across(starts_with("Trait"), ~ case_when(stri_replace_last_regex(., "_[0-9]", "") == "Size" ~ "Size\n(mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_curved_setaeP1P5_default" ~ "Sum of curved\nsetae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_spinesP6P7" ~ "Sum of spines",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Size" ~ "Standardised\nsize (mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_curved_setaeP1P5_default" ~ "Standardised sum\nof curved setae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_spinesP6P7" ~ "Standardised sum\nof spines",
                                                  TRUE ~ NA) %>% 
           factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines", "Standardised\nsize (mm)", "Standardised sum\nof curved setae", "Standardised sum\nof spines"))) ,
         Onset = case_when(Onset == "0" ~ "1\nMating success",
                           Onset == "1" ~ "2\nFecundity per mate",
                           Onset == "2" ~ "3\nOffspring survival") %>%
           factor(levels = c("1\nMating success", "2\nFecundity per mate",
                             "12\nFecundity", "3\nOffspring survival", "123\nReproductive success")))
         
var_covar_plot <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                                         traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines"))
var_covar_plot
ggsave(plot = var_covar_plot, "../Figures/tmp/Sexual_selection/Variance_covariance.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

var_covar_plot_std <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                                             traits = c("Standardised\nsize (mm)", "Standardised sum\nof curved setae", "Standardised sum\nof spines"))
var_covar_plot_std
ggsave(plot = var_covar_plot_std, "../Figures/tmp/Sexual_selection/Variance_covariance_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

########################################################
############ Plots specific cases gradients ############
########################################################
save_plot_regression(fitness_components, "Male", "praehirsuta", "Sum_spinesP6P7", "w2")

fitness_components %>% 
  filter(Sex == "Male", Species == "forsmani") %>% 
  plot_phenotype_fitness_regression("Size", "w123")
