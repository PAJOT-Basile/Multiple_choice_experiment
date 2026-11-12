########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "ggh4x", "parallel", "scales", "ggpubr", "stringi", "patchwork")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################################
################### Useful variables ###################
########################################################
set.seed(123456789)
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))

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

non_linear_differential <- function(var_b4, var_after, linear_differential){
  return(var_after - var_b4 + tcrossprod(linear_differential))
}

recompute_relative_fitnesses <- function(df){
  df$p0 <- 1 / nrow(df)
  df$w1 <- ifelse(is.na(df$W1 / mean(df$W1, na.rm = TRUE)), 0, df$W1 / mean(df$W1, na.rm = TRUE))
  df$p1 <- ifelse(is.na(df$p0 * df$w1), 0, df$p0 * df$w1)
  df$w2 <- ifelse(is.na(df$W2 / sum(df$p1 * df$W2, na.rm = TRUE)), 0, df$W2 / sum(df$p1 * df$W2, na.rm = TRUE))
  df$p2 <- ifelse(is.na(df$p1 * df$w2), 0, df$p1 * df$w2)
  df$w3 <- ifelse(is.na(df$W3 / sum(df$p2 * df$W3, na.rm = TRUE)), 0, df$W3 / sum(df$p2 * df$W3, na.rm = TRUE))
  df$p3 <- ifelse(is.na(df$p2 * df$w3), 0, df$p2 * df$w3)
  
  return(df)
}

compute_beta_and_S_in_bootstraps <- function(i, df, cols_to_take, onsets = c(1, 2, 3, 12, 123), resample = TRUE){
  if (resample){
    lines_to_sample <- sample(1:nrow(df), replace = TRUE, size = nrow(df))
  }else{
    lines_to_sample <- 1:nrow(df)
  }
  df_samp <- df[lines_to_sample, ]
  
  # Calculate the parameters
  df_samp <- recompute_relative_fitnesses(df_samp)
  
  
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
  C <- data.frame()
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
    
    p_1 <- df_samp %>% 
      select(paste0("p", onset_nb_1)) %>% 
      as.matrix()
    P_1 <- weighted_cov(Z, p_1)
      
    keep <- which(diag(P_1) > 0)
    
    if(nrow(Z) == 0 | length(keep) == 0){
      beta_vect[counter + 1:nb_traits] <- beta
    }else{
      Pk_1 <- P_1[keep, keep, drop = FALSE]
      
      ev <- eigen(Pk_1, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) < 1e-8 * max(ev)){
        beta_vect[counter + 1:nb_traits] <- beta   # relative, scale-free test
      } else {
        traits_to_keep <- colnames(Pk_1)
        
        res_beta <- tryCatch(as.vector(solve(Pk_1, S_vals[S_vals$Trait %in% traits_to_keep, corresponding_S_val])),
                        error = function(e) NULL)
        
        if (is.null(res_beta) || any(!is.finite(res_beta))){
          beta_vect[counter + 1:nb_traits] <- NA
        } else {
          beta_vect[counter + keep] <- res_beta
        }
        # Compute the non-linear selection differential
        p <- df_samp %>% 
          select(paste0("p", str_split_last(onset, ""))) %>% 
          as.matrix()
        P <- weighted_cov(Z, p)
        Pk <- P[rownames(P) %in% traits_to_keep, colnames(P) %in% traits_to_keep]
        
        C <- non_linear_differential(var_b4 = Pk_1,
                                     var_after = Pk,
                                     linear_differential = S_vals[S_vals$Trait %in% traits_to_keep, corresponding_S_val]) %>% 
            as.data.frame() %>% 
            rownames_to_column("Trait1") %>% 
            pivot_longer(!Trait1, names_to = "Trait2", values_to = "Stat_value") %>% 
            rowwise() %>% 
            mutate(comb_trait = paste(sort(c(Trait1, Trait2)), collapse = "_"),
                   Stat_name = paste0("C", onset)) %>%
            distinct(comb_trait, .keep_all = TRUE) %>% 
            select(-comb_trait) %>% 
          rbind(C)
      }
    counter <- counter + nb_traits
    }
  }
  names(beta_vect) <- paste0(names_beta_vect$Var1, "_", names_beta_vect$Var2)
  
  beta_vect %>%
    t() %>%
    as.data.frame() %>%
    mutate(Stat_name = "Beta") %>%
    rbind(S_vals %>%
            pivot_longer(-Trait, names_to = "Onset", values_to = "S") %>%
            mutate(Onset = paste(Trait, str_remove_all(Onset, "S"), sep = "_")) %>%
            select(-Trait) %>%
            pivot_wider(names_from = Onset, values_from = S) %>%
            mutate(Stat_name = "S")) %>%
    pivot_longer(!Stat_name, names_to = "Trait1", values_to = "Stat_value") %>%
    mutate(Trait2 = NA) %>% 
    rbind(C) %>% 
    return()
}

bootstap_S_and_grads <- function(sex, species, df, cols_to_take, nb_boots){
  nb_cores_to_use <- min(detectCores() - 1, nb_boots)
  cl <- makeCluster(nb_cores_to_use)
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, varlist = c("weighted_cov", "recompute_relative_fitnesses", "non_linear_differential", "str_split_last"))
  
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

compute_covar_matrices <- function(onset_nb, cols_to_take, species, sex, phenos){
  onset_nb_1 <- case_when(onset_nb %in% c(1, 2, 3, 4) ~ onset_nb - 1,
                          onset_nb %in% c(12, 123) ~ 0)
  Z <- phenos %>% 
    select(all_of(cols_to_take))  %>%
    select_if(~ var(., na.rm = TRUE) != 0) %>% 
    as.matrix()
  p <- phenos %>% 
    select(paste0("p", onset_nb_1)) %>% 
    as.matrix()
  
  P <- weighted_cov(Z, p) 
  keep <- which(diag(P) > 0)

  if(nrow(Z) == 0 | length(keep) == 0){
    P_i <- data.frame(
      "Trait1" = NA,
      "Trait2" = NA,
      "var_covar" = NA,
      "Onset" = onset_nb_1
    )
  }else{
    Pk <- P[keep, keep, drop = FALSE]
    if (length(Pk) == 1){
      P_i <- Pk %>% 
        as_tibble() %>% 
        mutate(Trait1 = names(which(colSums(P) != 0)),
               Trait2 = names(which(colSums(P) != 0)),
               Onset = onset_nb_1) %>% 
        rename(var_covar = !!sym(names(which(colSums(P) != 0))))
      }else{
        P_i <- Pk %>% 
          as.data.frame() %>% 
          rownames_to_column("Trait1") %>% 
          pivot_longer(!Trait1, names_to = "Trait2", values_to = "var_covar") %>% 
          mutate(Onset = onset_nb_1)
      }
  }
  return(P_i)
}

plot_linerange_traits <- function(df, traits, stat_name, ylab = ""){
  
  bold_labels <- ifelse(levels(df$Onset) %in% c("12\nNumber of offspring", "123\nReproductive success"), "bold", "plain")
  x_text_colour <- ifelse(levels(df$Onset) %in% c("12\nNumber of offspring", "123\nReproductive success"), "black", "grey30")
  
  forsmani <- df %>%
    filter(Species == "J. forsmani",
           Trait1 %in% traits,
           grepl(stat_name, Stat_name),
           !(Sex == "Female" & grepl("Sum", Trait1) & grepl(stat_name, Stat_name)),
           !is.na(Stat_value),
           ) 
  
  if (stat_name == "C"){
    forsmani <- forsmani %>% 
      filter(Trait2 %in% traits)
  }
  
  forsmani <- forsmani %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), colour = "black", lwd = 1.2) +
    geom_linerange(aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5)) +
    geom_linerange(data = . %>% 
                     filter(Onset %in% c("12\nNumber of offspring", "123\nReproductive success")),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1.4) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    geom_text(aes(x = Onset, y = Inf, label = Signif), size = 8, vjust = 1) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) 
  
  
  if (stat_name == "C"){
    forsmani <- forsmani + 
      facet_grid2(Sex + Trait1 + Trait2 ~ ., scales = "free", independent = "y", strip = strip_nested(
        background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
        by_layer_y = TRUE)) +
      facetted_pos_scales(y = list(
        Trait1 == "Size\n(mm)" & stat_name == "Beta" ~ scale_y_continuous(labels = function(x){sprintf("%.1f", x)}),
        TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
      ))
  }else{
    forsmani <- forsmani +
      facet_grid2(Sex + Trait1 ~ ., scales = "free", independent = "y", strip = strip_nested(
        background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
        by_layer_y = TRUE)) +
    facetted_pos_scales(y = list(
      Trait1 == "Size\n(mm)" & stat_name == "Beta" ~ scale_y_continuous(labels = function(x){sprintf("%.1f", x)}),
      TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
    ))
  }
  forsmani <- forsmani +
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
           Trait1 %in% traits,
           grepl(stat_name, Stat_name),
           !(Sex == "Female" & grepl("Sum", Trait1) & grepl(stat_name, Stat_name)),
           !is.na(Stat_value),
    ) 
  
  if (stat_name == "C"){
    praehirsuta <- praehirsuta %>% 
      filter(Trait2 %in% traits)
  }
  
  praehirsuta <- praehirsuta %>%
    ggplot() +
    geom_hline(aes(yintercept = 0), colour = "black", lwd = 1.2) +
    geom_linerange(aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5)) +
    geom_linerange(data = . %>% 
                     filter(Onset %in% c("12\nNumber of offspring", "123\nReproductive success")),
                   aes(x = Onset, ymin = CI_2.5, ymax = CI_97.5), lwd = 1.4) +
    geom_point(aes(x = Onset, y = Stat_value, fill = Species), size = 4,
               show.legend = TRUE, colour = "black", pch = 21) +
    geom_text(aes(x = Onset, y = Inf, label = Signif), size = 8, vjust = 1) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      drop = FALSE,
                      breaks = c("J. forsmani", "J. praehirsuta")) 
  
  
  if (stat_name == "C"){
    praehirsuta <- praehirsuta + 
      facet_grid2(Sex + Trait1 + Trait2 ~ ., scales = "free", independent = "y", strip = strip_nested(
        background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
        by_layer_y = TRUE)) +
      facetted_pos_scales(y = list(
        Trait1 == "Size\n(mm)" & stat_name == "Beta" ~ scale_y_continuous(labels = function(x){sprintf("%.1f", x)}),
        TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
      ))
  }else{
    praehirsuta <- praehirsuta +
      facet_grid2(Sex + Trait1 ~ ., scales = "free", independent = "y", strip = strip_nested(
        background_y = list(element_rect(fill = "white"), element_rect(fill = "grey85")),
        by_layer_y = TRUE)) +
      facetted_pos_scales(y = list(
        Trait1 == "Size\n(mm)" & stat_name == "Beta" ~ scale_y_continuous(labels = function(x){sprintf("%.1f", x)}),
        TRUE ~ scale_y_continuous(labels = function(x) sprintf("%.2f", x))
      ))
  }
  praehirsuta <- praehirsuta +
    labs(y = ylab) +
    theme_grid  +
    theme(legend.text = element_text(face = "italic"),
          strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(face = bold_labels,
                                     colour = x_text_colour))
  

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
  fitness_onset <- c("b4" = ifelse(fitness %in% c("w12", "w123", "W12", "W123"), 0, as.numeric(str_split_fixed(fitness, "", 2)[, 2]) - 1),
                     "after" = as.numeric(str_split_last(fitness, "")))

  phenotype <- ifelse(grepl("setae", phenotype), paste0(phenotype, "_default"), phenotype)
  
  phenotypes <- c("b4" = paste0(phenotype, "_", fitness_onset["b4"]),
                  "after" = paste0(phenotype, "_", fitness_onset["after"]))
  
  ps <- c("b4" = unname(ifelse(fitness_onset["b4"] == 0, paste0("p", fitness_onset["b4"]), paste0("p", fitness_onset["b4"] - 1))),
          "after" = paste0("p", as.numeric(str_split_last(fitness, "")) - 1))
  
  
  relative_fitness <- c("b4" = unname(ifelse(fitness_onset["b4"] == 0, NA, paste0("w", fitness_onset["b4"]))),
                        "after" = paste0("w", str_split_last(fitness, "")))

  df <- recompute_relative_fitnesses(df)

  if (fitness_onset["b4"] == 0){
    df[[phenotypes["b4"]]] <- df[[phenotype]] * df[[ps["b4"]]]

  }else{
    df <- df %>% 
      filter(!!sym(relative_fitness["b4"]) != 0)
    df[[phenotypes["b4"]]] <- df[[phenotype]] * df[[ps["b4"]]] * df[[relative_fitness["b4"]]]
  }

  df[[phenotypes["after"]]] <- df[[phenotype]] * df[[ps["after"]]] * df[[relative_fitness["after"]]]
  
  start_arrow <- sum(df[[phenotypes["b4"]]], na.rm = TRUE)
  end_arrow <- sum(df[[phenotypes["after"]]], na.rm = TRUE)
  height_arrow <- 1.5*mean(df[[fitness]], na.rm = TRUE)
  
  sum_lm <- df %>% 
    filter(!is.na(!!sym(phenotype)), !is.na(!!sym(fitness))) %>% 
    lm(as.formula(paste(fitness, "~", phenotype)), data = .) %>% 
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
    filter(!is.na(!!sym(phenotype)), !is.na(!!sym(fitness))) %>% 
    ggplot() +
    geom_point(aes(x = !!sym(phenotype), y = !!sym(fitness)), colour = "black", size = 4, alpha = 0.5) +
    geom_vline(aes(xintercept = sum(!!sym(phenotypes["b4"]), na.rm = TRUE)),
               colour = "black", lwd = 1.1) +
    geom_vline(aes(xintercept = sum(!!sym(phenotypes["after"]), na.rm = TRUE)),
               colour = "grey40", lwd = 1.1) +
    geom_smooth(aes(x = !!sym(phenotype), y = !!sym(fitness)), method = "lm",
                se = FALSE, colour = "grey40", lty = "dashed") +
    geom_smooth(aes(x = !!sym(phenotype), y = !!sym(fitness)), method = "lm",
                se = FALSE, colour = "grey40", formula = y ~ poly(x, degree = 2),
                lty = "dotted") +
    annotate("segment", x = start_arrow, xend = end_arrow, y = height_arrow,
             arrow = arrow(length = unit(0.5, "cm")), lwd = 1.1) +
    # annotate("text", x = Inf, y = Inf, label = lab,
    #          hjust = 1.2, vjust = 1.5, size = 7, parse = TRUE) +
    scale_y_continuous(limits = c(0, NA)) +
    my_theme
  return(p)
}

save_plot_regression <- function(df, sex, species, pheno, fitness){
  y_text <- case_when(fitness %in% c("w1", "W1") ~ "1\nMating success",
                      fitness %in% c("w2", "W2") ~ "2\nNumber of offspring\nper mate",
                      fitness %in% c("w12", "W12") ~ "12\nNumber of offspring",
                      fitness %in% c("w3", "W3") ~ "3\nOffspring survival",
                      fitness %in% c("w123", "W123") ~ "123\nReproductive success") %>% 
    factor(levels = c("1\nMating success", "2\nNumber of offspring\nper mate",
                      "12\nNumber of offspring", "3\nOffspring survival", "123\nReproductive success"))
  x_text <- case_when(pheno == "Size" ~ "Size\n(mm)",
                      pheno == "Sum_curved_setaeP1P5" ~ "Sum of curved\nsetae",
                      pheno == "Sum_spinesP6P7" ~ "Sum of spines\n",
                      pheno == "St_Size" ~ "Standardised size\n(mm)",
                      pheno == "St_Sum_curved_setaeP1P5" ~ "Standardised sum of curved\nsetae",
                      pheno == "St_Sum_spinesP6P7" ~ "Standardised sum of spines") %>% 
    factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n", "Standardised size\n(mm)",
                      "Standardised sum of curved\nsetae", "Standardised sum of spines"))
  
  p <- df %>% 
    filter(Sex == sex, Species == species) %>% 
    plot_phenotype_fitness_regression(pheno, fitness) +
    geom_hline(yintercept = 0, colour = "black", lwd = 1) +
    labs(
      x = x_text,
      y = y_text,
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
  left_join(small_metadata,
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


traits_to_use <- matrix(c("Size", "Sum_curved_setaeP1P5_default", "Sum_spinesP6P7",
                          "St_Size", "St_Sum_curved_setaeP1P5_default", "St_Sum_spinesP6P7"),
                        nrow = 2, ncol = 3, byrow = TRUE)


selection_differentials_gradients <- lapply(c("praehirsuta", "forsmani"), function(species, df, traits_to_use){
  df_filt <- df %>% 
    filter(Species == species)
  
  lapply(c("Male", "Female"), function(sex, species, df, traits_to_use){
    df_filt <- df %>% 
      filter(Sex == sex)

    apply(traits_to_use, 1, function(cols_to_take, sex, species, df){

      compute_beta_and_S_in_bootstraps(1, df, cols_to_take, resample = FALSE)

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
  mutate(Onset = case_when(Stat_name %in% c("Beta", "S") ~ str_split_last(Trait1, "_"),
                           TRUE ~ str_remove_all(Stat_name, "C")),
         Trait1 = case_when(Stat_name %in% c("Beta", "S") ~ str_remove(Trait1, "_[^_]*$"),
                            TRUE ~ Trait1)) %>% 
  mutate(Stat_name = case_when(Stat_name %in% c("Beta", "S") ~ paste(Stat_name, Onset, sep = "_"),
                               TRUE ~ str_replace_all(Stat_name, "(?<=[A-Z])(?=\\d)", "_")))

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
#       bootstap_S_and_grads(sex, species, df, cols_to_take, nb_boots)
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
#   mutate(Stat_name = case_when(Stat_name %in% c("Beta", "S") ~ paste(Stat_name, str_split_last(Trait1, "_"), sep = "_"),
#                                TRUE ~ str_replace_all(Stat_name, "(?<=[A-Z])(?=\\d)", "_")),
#          Trait1 = case_when(grepl("Beta|S", Stat_name) ~ str_remove(Trait1, "_[^_]*$"),
#                             TRUE ~ Trait1))
# 
# 
# distrib_sel_diffs_and_grads <- distrib_sel_diffs_and_grads %>% 
#   ungroup()
# save(distrib_sel_diffs_and_grads, file = "Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000_three_episodes.rdata")
load("Associated_data/Rdata/distrib_selection_differentials_and_gradients_bootstraps_50000_three_episodes.rdata")


CI_stats <- distrib_sel_diffs_and_grads %>% 
  group_by(Species, Sex, Trait1, Trait2, Stat_name) %>% 
  summarize(CI_2.5 = quantile(Stat_value, probs = 0.025, na.rm = TRUE),
            CI_97.5 = quantile(Stat_value, probs = 0.975, na.rm = TRUE),
            mean_val_boots = mean(Stat_value, na.rm = TRUE)) %>% 
  ungroup()

signif <- distrib_sel_diffs_and_grads %>%
  group_by(Sex, Species, Trait1, Trait2, Stat_name) %>%
  proba_different_from_0_with_distribution(Stat_value) %>% 
  ungroup() %>% 
  mutate(Signif = case_when(pval < 0.05 & pval >= 0.01 ~ "*",
                            pval < 0.01 & pval > 0.001 ~ "**",
                            pval <= 0.001 ~ "***",
                            TRUE ~ ""))

rm(distrib_sel_diffs_and_grads)
########################################################
##################### Plot results #####################
########################################################
selection_differentials_gradients <- selection_differentials_gradients %>%
  left_join(signif, by = c("Sex", "Species", "Trait1", "Trait2", "Stat_name")) %>%
  left_join(CI_stats, by = c("Sex", "Species", "Trait1", "Trait2", "Stat_name")) %>% 
  mutate(Onset = case_when(Onset == "1" ~ "1\nMating success",
                           Onset == "2" ~ "2\nNumber of offspring\nper mate",
                           Onset == "12" ~ "12\nNumber of offspring",
                           Onset == "3" ~ "3\nOffspring survival",
                           Onset == "123" ~ "123\nReproductive success") %>% 
           factor(levels = c("1\nMating success", "2\nNumber of offspring\nper mate",
                             "12\nNumber of offspring", "3\nOffspring survival", "123\nReproductive success")),
         across(starts_with("Trait"), ~ case_when(. == "Size" ~ "Size\n(mm)",
                           . == "Sum_curved_setaeP1P5_default" ~ "Sum of curved\nsetae",
                           . == "Sum_spinesP6P7" ~ "Sum of spines\n",
                           . == "St_Size" ~ "Standardised size\n(mm)",
                           . == "St_Sum_curved_setaeP1P5_default" ~ "Standardised sum of curved\nsetae",
                           . == "St_Sum_spinesP4P7" ~ "Standardised sum of spines") %>% 
           factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n",
                             "Standardised size\n(mm)", "Standardised sum of curved\nsetae",
                             "Standardised sum of spines"))),
         Species = paste0("J. ", Species) %>% factor(levels = c("J. forsmani", "J. praehirsuta")))


# Look at the selection differentials
selection_differentials_plot <- plot_linerange_traits(df = selection_differentials_gradients,
                                                      traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n"),
                                                      stat_name = "S",
                                                      ylab = "Selection differentials")

selection_differentials_plot
ggsave(plot = selection_differentials_plot, "../Figures/Suppfigs/Selection_differentials.png",
       grDevices::png, scale = 3, width = 1200, height = 1300, units = "px")


selection_gradients_plot <- plot_linerange_traits(df = selection_differentials_gradients,
                                                  traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n"),
                                                  stat_name = "Beta",
                                                  ylab = "Selection gradients")
selection_gradients_plot
ggsave(plot = selection_gradients_plot, "../Figures/Suppfigs/Selection_gradients.png",
       grDevices::png, scale = 3, width = 1200, height = 1300, units = "px")



non_linear_seletcion_differentials <- plot_linerange_traits(df = selection_differentials_gradients,
                                                            traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n"),
                                                            stat_name = "C",
                                                            ylab = "Non-linear selection differentials")


non_linear_seletcion_differentials
ggsave(plot = non_linear_seletcion_differentials, "../Figures/Suppfigs/Non_linear_selection_differentials.png",
       grDevices::png, scale = 3, width = 1200, height = 1800, units = "px")


# Look at the standardised traits
selection_differentials_plot_std <- plot_linerange_traits(df = selection_differentials_gradients,
                                                          traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                          stat_name = "S",
                                                          ylab = "Selection differentials")
selection_differentials_plot_std
ggsave(plot = selection_differentials_plot_std, "../Figures/tmp/Standardised_traits/Selection_differentials_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

selection_gradients_plot_std <- plot_linerange_traits(df = selection_differentials_gradients,
                                                      traits = c("Standardised size\n(mm)", "Standardised sum of curved\nsetae", "Standardised sum of spines"),
                                                      stat_name = "Beta",
                                                      ylab = "Selection gradients")
selection_gradients_plot_std
ggsave(plot = selection_gradients_plot_std, "../Figures/tmp/Standardised_traits/Selection_gradients_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")


########################################################
################## Covariance matrices ##################
########################################################

variance_covariance_matrices <- apply(traits_to_use, 1, function(cols_to_take, df){
  lapply(c("forsmani", "praehirsuta"), function(species, df, cols_to_take){
    
    df_filt <- df %>% 
      filter(Species == species)
    lapply(c("Male", "Female"), function(sex, species, df, cols_to_take){
      
      df_filt <- df %>% 
        filter(Sex == sex)
      
      lapply(c(1, 2, 3), compute_covar_matrices, cols_to_take, species, sex, df_filt) %>% 
        bind_rows() %>% 
        mutate(Sex = sex,
               Species = species)
      
    }, species, df_filt, cols_to_take) %>% 
      bind_rows()
    
  }, df, cols_to_take) %>% 
    bind_rows()
}, fitness_components %>%
  select(-contains("trunc"))) %>% 
  bind_rows() %>% 
  mutate(across(starts_with("Trait"), ~ case_when(stri_replace_last_regex(., "_[0-9]", "") == "Size" ~ "Size\n(mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_curved_setaeP1P5_default" ~ "Sum of curved\nsetae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "Sum_spinesP6P7" ~ "Sum of spines\n",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Size" ~ "Standardised\nsize (mm)",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_curved_setaeP1P5_default" ~ "Standardised sum\nof curved setae",
                                                  stri_replace_last_regex(., "_[0-9]", "") == "St_Sum_spinesP6P7" ~ "Standardised sum\nof spines",
                                                  TRUE ~ NA) %>% 
           factor(levels = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n", "Standardised\nsize (mm)", "Standardised sum\nof curved setae", "Standardised sum\nof spines"))) ,
         Onset = case_when(Onset == "0" ~ "1\nMating success",
                           Onset == "1" ~ "2\nNumber of offspring\nper mate",
                           Onset == "2" ~ "3\nOffspring survival") %>%
           factor(levels = c("1\nMating success", "2\nNumber of offspring\nper mate",
                             "12\nNumber of offspring", "3\nOffspring survival")))
         
var_covar_plot <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                                         traits = c("Size\n(mm)", "Sum of curved\nsetae", "Sum of spines\n"))
var_covar_plot +
  scale_fill_gradientn(trans = "log",
                       colours = c("navy", "dodgerblue", "skyblue", "white", "orange2", "red", "firebrick"))
ggsave(plot = var_covar_plot, "../Figures/Suppfigs/Variance_covariance.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

var_covar_plot_std <- plot_heatmap_var_covar(df = variance_covariance_matrices,
                                             traits = c("Standardised\nsize (mm)", "Standardised sum\nof curved setae", "Standardised sum\nof spines"))
var_covar_plot_std
ggsave(plot = var_covar_plot_std, "../Figures/tmp/Standardised_traits/Variance_covariance_standardised_traits.png",
       grDevices::png, scale = 4, width = 1200, height = 1300, units = "px")

########################################################
############ Plots specific cases gradients ############
########################################################
# Trait nmaes:
#   Size
#   Sum_curved_setaeP1P5
#   Sum_spinesP6P7
fitness_components %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  ggplot(aes(x = Size, y = W1, colour = Species)) +
  geom_point(colour = "black", alpha = 0.5, size = 3) +
  geom_smooth(colour = "grey40", method = "lm", lty = "dashed", se = FALSE) +
  geom_smooth(colour = "blue", method = "lm", formula = y ~ poly(x, degree = 2), se = FALSE) +
  my_theme


fitness_components %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  lm(data = ., w1 ~ Size + I(Size^2)) %>%
  # lm(data = ., w1 ~ Size) %>% 
  AIC()


fitness_components %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  plot_phenotype_fitness_regression("Size", "w2")


fitness_components %>% 
  filter(Sex == "Male", Species == "praehirsuta") %>% 
  ggplot() + 
  geom_point(aes(x = Sum_curved_setaeP1P5_default, y = w3, colour = Sum_spinesP6P7),
             size = 4) +
  my_theme


save_plot_regression(fitness_components, "Male", "praehirsuta", "Sum_spinesP6P7", "w3")
  


##############################################################################
plot_figure_trait_fitness_differentials <- function(s_values, fitness_components, species, sex, trait, tag_name = ""){
  xmas <- s_values %>%
    filter(Species == species, Sex == sex,
           Trait1 == trait,
           grepl("S", Stat_name)) %>% 
    pivot_longer(cols = contains("CI_"), names_to = "CI", values_to = "CI_val") %>% 
    pull(CI_val) %>% 
    abs() %>% 
    max(na.rm = TRUE)
  
  xmac <- s_values %>%
    filter(Species == species, Sex == sex,
           Trait1 == trait,
           Trait2 == trait,
           grepl("C", Stat_name)) %>% 
    pivot_longer(cols = contains("CI_"), names_to = "CI", values_to = "CI_val") %>% 
    pull(CI_val) %>% 
    abs() %>% 
    max(na.rm = TRUE)
  
  selection_differentials_s <- s_values %>%
    mutate(Onset = factor(Onset, levels = c("123\nReproductive success", "3\nOffspring survival",
                                            "12\nNumber of offspring", "2\nNumber of offspring\nper mate",
                                            "1\nMating success"))) %>% 
    filter(Species == species, Sex == sex,
           Trait1 == trait,
           grepl("S", Stat_name)) %>% 
    ggplot() +
    geom_vline(aes(xintercept = 0), colour = "black", lwd = 0.9, lty = "dashed") +
    geom_linerange(aes(xmin = CI_2.5, xmax = CI_97.5, y = Onset)) +
    geom_point(aes(x = Stat_value, y = Onset, fill = Species), size = 5,
               show.legend = TRUE, colour = "black", pch = 21) +
    geom_text(aes(x = ifelse(abs(CI_97.5) > abs(CI_2.5), CI_97.5, CI_2.5),
                  y = Onset, label = Signif), size = 10, vjust = -0.5) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      guide = "none") +
    labs(x = "Linear\nselection differentials") +
    xlim(-xmas, xmas) +
    my_theme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.ticks.length.x = unit(1, "mm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  selection_differentials_c <- s_values %>%
    mutate(Onset = factor(Onset, levels = c("123\nReproductive success", "3\nOffspring survival",
                                            "12\nNumber of offspring", "2\nNumber of offspring\nper mate",
                                            "1\nMating success"))) %>% 
    filter(Species == species, Sex == sex,
           Trait1 == trait,
           Trait2 == trait,
           grepl("C", Stat_name)) %>% 
    ggplot() +
    geom_vline(aes(xintercept = 0), colour = "black", lwd = 0.9, lty = "dashed") +
    geom_linerange(aes(xmin = CI_2.5, xmax = CI_97.5, y = Onset)) +
    geom_point(aes(x = Stat_value, y = Onset, fill = Species), size = 5,
               show.legend = TRUE, colour = "black", pch = 21) +
    geom_text(aes(x = ifelse(abs(CI_97.5) > abs(CI_2.5), CI_97.5, CI_2.5),
                  y = Onset, label = Signif), size = 10, vjust = -0.5) +
    scale_fill_manual(values = c("J. praehirsuta" = "navy", "J. forsmani" = "#3A9AB2"),
                      guide = "none") +
    labs(x = "Non-linear\nselection differentials") +
    xlim(-xmac, xmac) +
    my_theme +
    theme(legend.text = element_text(face = "italic"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.ticks.length.x = unit(1, "mm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  sp <- str_remove_all(species, "J. ")
  tr <- case_when(trait == "Size\n(mm)" ~ "Size",
                  trait == "Sum of curved\nsetae" ~ "Sum_curved_setaeP1P5",
                  trait == "Sum of spines\n" ~ "Sum_spinesP6P7")
  
  tra <- ifelse(tr == "Sum_curved_setaeP1P5", "Sum_curved_setaeP1P5_default", tr)
  
  xmin_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp) %>% 
    pull(!!sym(tra)) %>% 
    min(na.rm = TRUE)
  
  xmax_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp) %>% 
    pull(!!sym(tra)) %>% 
    max(na.rm = TRUE)
    
  
  w1_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp) %>% 
    plot_phenotype_fitness_regression(tr, "w1") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(xmin_trait, xmax_trait)) +
    labs(y = "w1\nRelative mating success",
         tag = tag_name) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.tag = element_text(size = 40))
  
  w2_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp,
           W1 != 0) %>% 
    plot_phenotype_fitness_regression(tr, "w2") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(xmin_trait, xmax_trait)) +
    labs(y = "w2\nRelative number of offspring per mate") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  w12_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp) %>% 
    plot_phenotype_fitness_regression(tr, "w12") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(xmin_trait, xmax_trait)) +
    labs(y = "w12\nRelative number of offspring") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  w3_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp,
           W2 != 0) %>% 
    plot_phenotype_fitness_regression(tr, "w3") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(xmin_trait, xmax_trait)) +
    labs(y = "w3\nRelative offspring survival") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  w123_trait <- fitness_components %>%
    filter(Sex == sex, Species == sp) %>% 
    plot_phenotype_fitness_regression(tr, "w123") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(xmin_trait, xmax_trait)) +
    labs(y = "w123\nRelative reproductive success",
         x = trait) +
    theme(axis.text.x = element_text(size = 18),
          axis.title.x = element_text(size = 20))
  
  text_annotation <- paste0(sex, "~italic(", str_replace_all(species, " ", "~"), ")~on~",
                            str_remove_all(trait, "\n$") %>% 
                              str_replace_all(., "\n", "~") %>% 
                              str_replace_all(., " ", "~"))
  
  ann1 <-  ggplot() +
    geom_text(aes(x = 0, y = 0), label = text_annotation,
              parse = TRUE, size = 8, vjust = 1) +
    theme_void() 
    
  
  
  trait_fitness_plot <- ann1 +
    w1_trait + selection_differentials_s + selection_differentials_c +
    w2_trait + 
    w12_trait +
    w3_trait + 
    w123_trait +
    plot_layout(design = "
              AAA
              BCD
              ECD
              FCD
              GCD
              HCD",
                widths = c(0.85, 0.25, 0.25),
                heights = c(0.15, 1, 1, 1, 1, 1))
  
  return(trait_fitness_plot)
  
}


#' Traits:
#' "Size\n(mm)"
#' "Sum of spines\n"
#' "Sum of curved\nsetae"
plot_figure_trait_fitness_differentials(s_values = selection_differentials_gradients,
                                        fitness_components = fitness_components,
                                        species = "J. praehirsuta",
                                        sex = "Male",
                                        trait = "Sum of curved\nsetae")

size_male_forsmani <- plot_figure_trait_fitness_differentials(s_values = selection_differentials_gradients,
                                                              fitness_components = fitness_components,
                                                              species = "J. forsmani",
                                                              sex = "Male",
                                                              trait = "Size\n(mm)",
                                                              tag_name = "(A)")

spines_male_forsmani <- plot_figure_trait_fitness_differentials(selection_differentials_gradients, fitness_components,
                                                                "J. forsmani", "Male", "Sum of spines\n",
                                                                tag_name = "(B)")

forsmani_plot_s <- ggarrange(size_male_forsmani, spines_male_forsmani,
                             ncol = 2)

ggsave(forsmani_plot_s, 
       filename = "../Figures/Selection_differentials_with_plots_forsmani.png",
       scale = 4, height = 1800, width = 2600, units = "px")



size_female_prae <- plot_figure_trait_fitness_differentials(selection_differentials_gradients, fitness_components,
                                                                "J. praehirsuta", "Female", "Size\n(mm)",
                                                            tag_name = "(A)")
setae_male_prae <- plot_figure_trait_fitness_differentials(selection_differentials_gradients, fitness_components,
                                                           "J. praehirsuta", "Male", "Sum of curved\nsetae",
                                                           tag_name = "(B)")

praehirsuta_plot_s <- ggarrange(size_female_prae, setae_male_prae,
                             ncol = 2)

ggsave(praehirsuta_plot_s, 
       filename = "../Figures/Selection_differentials_with_plots_praehirsuta.png",
       scale = 4, height = 1800, width = 2600, units = "px")


######################## Test Multivariate normal distribution of traits ########################
x <- metadata %>%
  filter(Sex == "M", Species == "forsmani") %>%
  select(Sum_curved_setaeP1P5_default, Sum_spinesP6P7, Size)

 
plotly::plot_ly(x,
                x = ~ Size,
                y = ~ Sum_curved_setaeP1P5_default,
                z = ~ Sum_spinesP6P7,
                xlab = "Size")
