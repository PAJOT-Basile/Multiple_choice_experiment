#' Script designed to run the animal model to compute the heritability of sexual
#' traits (size of individuals and for males, the number of sexual setae and
#' number of spines)
#' 
#' Basile Pajot
#' Last updated: 17/11/2026

########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "MCMCglmm", "parallel", "patchwork")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################
########### Useful variables ###########
########################################
set.seed(123456789)
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20))

theme_grid <- my_theme +
  theme(plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
        legend.margin = margin(t = 15, r = 15, b = 15, l = 15),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 15),
        text = element_text(size = 12))

########################################
########### Useful functions ###########
########################################
source("98_Transversal_scripts/Transversal_functions.r")

run_model_on_modelid <- function(model_id, effect_traits, trait, sp, metadata, pedigree, compute_heritability = FALSE){
  print(paste0("       - Model_ID: ", model_id))
  
  fixed_effect_traits <- effect_traits %>%
    filter(Model_ID == model_id) %>%
    pull(Fixed_effects) %>% 
    unlist()
  random_effect_traits <- effect_traits %>%
    filter(Model_ID == model_id) %>%
    pull(Random_effects) %>% 
    unlist()
  # Define the family of the link function depending on the data input (the total
  # number of sexual setae or the total number of spines are taken in a poisson
  # distribution while the size of individuals is taken in a gaussian distribution).
  #
  # In addition, if the focal trait is either the number of spines or the number
  # of sexual setae, we remove the females from the data frame to not increase the
  # variance artificially.
  if (grepl("Sum", trait)){
    family_model <- "poisson"
    
    metadata_sp <- metadata %>%
      filter(Sex == "M")
  } else {
    family_model <- "gaussian"
    metadata_sp <- metadata
  }
  
  
  # Filter the metadata to remove the individuals that are not needed
  metadata_sp <- metadata_sp  %>%
    filter(animal %in% pedigree$animal,
           Species == sp)
  
  # We can also start to build the formulas that will be used in the model
  fixed_effect_formula <- paste(trait, "~", 1)
  random_effect_formula <- "~ animal"
  
  
  if (any(!is.na(fixed_effect_traits))){
    # Here, we remove the missing values contained in the fixed effect variables
    for (fix_effect in fixed_effect_traits){
      if (is.na(fix_effect)) next
      metadata_sp <- metadata_sp %>%
        filter(!is.na(!!sym(fix_effect)))
    }
    
    # Here, we add the other fixed effects in the formula
    fixed_effect_formula <- paste(fixed_effect_formula,  "+", paste(fixed_effect_traits, collapse = " + "))
  }
  
  # We do the same thing for the random effects
  if (any(!is.na(random_effect_traits))){
    # Here, we remove the missing values contained in the random effect variables
    for (rand_effect in random_effect_traits){
      if (is.na(rand_effect)) next
      metadata_sp <- metadata_sp %>%
        filter(!is.na(!!sym(rand_effect)))
    }
    # Here, we add the other random effects in the formula
    random_effect_formula <- paste(random_effect_formula,  "+", paste(random_effect_traits, collapse = " + "))
  }
  
  # Convert the filtered table to a data frame for the model to work
  metadata_sp <- metadata_sp %>%
    filter(!is.na(!!sym(trait))) %>%
    as.data.frame()
  
  
  # Compute the variances of all traits (focal trait as well as the traits in the
  # fixed and random effects)
  variances_traits <- lapply(as.vector(na.omit(c(trait, unlist(random_effect_traits)))), function(x, df){
    if (is.na(x) || x == "Sex") return(NULL)
    var(df[[x]]) %>%
      as.data.frame() %>%
      rename(Variance = ".") %>%
      mutate(Trait = x) %>%
      return()
  }, metadata_sp) %>%
    bind_rows() %>%
    mutate(Effect = case_when(Trait == trait ~ "Focal",
                              Trait %in% random_effect_traits ~ "Random",
                              TRUE ~ NA))
  
  # If the variance is larger than 1, we decide to use the default prior (V = 1, nu = 0.002),
  # otherwise, we will use the prior for the small variances: V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000
  priors_random <- lapply(as.vector(na.omit(c(trait, unlist(random_effect_traits)))), function(x, df){
    if (is.na(x) || x == "Sex"){
      prior_val <- list(V = 1, nu = 0.002)
    }else{
      v <- df %>%
        filter(Trait == x) %>%
        pull(Variance)
      if (v >= 1){
        prior_val <- list(V = 1, nu = 0.002)
      } else {
        prior_val <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)
      }
    }
    return(prior_val)
  }, variances_traits)
  
  names(priors_random) <- paste0("G", 1:length(as.vector(na.omit(c(trait, unlist(random_effect_traits))))))
  
  
  prior <- list(R = list(V = 1, nu = 0.002),
                G = priors_random)
  
  # Run the model
  model <- MCMCglmm(formula(fixed_effect_formula),
                    random = formula(random_effect_formula),
                    family = family_model,
                    prior = prior,
                    pedigree = pedigree,
                    data = metadata_sp,
                    verbose = FALSE,
                    nitt = 2000000,
                    burnin = 200000,
                    thin = 2000)
                    # nitt = 100,
                    # burnin = 10,
                    # thin = 10)
  
  
  if (compute_heritability){
    saveRDS(model,
            paste0("Associated_data/Heritability/data/Model_", trait, "_", sp, "_", model_id, ".rds"))
    
    herit <- model$VCV[, "animal"] / rowSums(model$VCV)
    effective_size_heritability <- effectiveSize(herit)
    h2 <- mean(herit)
  }
  
  
  # Get some output stats
  DIC_mod <- summary(model)$DIC
  effective_size_fixed_effects <- effectiveSize(model[["Sol"]])
  effective_size_random_effects <- effectiveSize(model[["VCV"]])[names(effectiveSize(model[["VCV"]])) != "units"]
  
  # Store the output in a data frame
  out_df <- data.frame("Model_ID" = model_id,
                       "Trait" = trait,
                       "Species" = sp,
                       "DIC" = DIC_mod,
                       "Effective_size_fixed_effects" = paste(names(effective_size_fixed_effects) %>% 
                                                                str_remove_all("\\(") %>%
                                                                str_remove_all("\\)"), effective_size_fixed_effects,
                                                              collapse = "; "),
                       "Effective_size_random_effects" = paste(names(effective_size_random_effects) %>% 
                                                                 str_replace_all("animal", "pedigree"), effective_size_random_effects,
                                                               collapse = "; "),
                       "Random_effect" = random_effect_formula,
                       "Fixed_effects" = fixed_effect_formula,
                       "Variance_trait" = variances_traits %>%
                         filter(Effect == "Focal") %>%
                         pull(Variance),
                       "Variance_random_effects" = paste(random_effect_traits, variances_traits %>%
                                                           filter(Effect == "Random") %>%
                                                           pull(Variance),
                                                         collapse = "; "))
  
  if (compute_heritability){
    out_df <- out_df %>% 
      mutate(Effective_size_heriatbility = effective_size_heritability,
             h2 = h2)
  }
  
  return(out_df)
}

########################################
############# Import data ##############
########################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b"),
         animal = as.integer(str_remove_all(ID_DNA_RAD, "E"))) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(  # For the offspring, the number of setae was counted more precisely than for the parents, 
    # especially on the carps of the different pereiopods. For the values to be 
    # comparable between parents and offspring, we remove the carps from calculations.
    Sum_curved_setaeP1P5_trunc = Sum_curved_setaeP1P5_default - default_P1.carp.curv.setae - default_P2.carp.curv.setae - default_P3.carp.curv.setae - default_P4.carp.curv.setae - default_P5.carp.curv.setae,
    across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup() %>% 
  filter(!is.na(Sex))

# Make the pedigree
load("Associated_data/Rdata/Families_first_brood.rdata")
pedigree <- Families_first_brood %>%
  select(Mom) %>%
  rename(animal = Mom) %>%
  unique() %>%
  mutate(sire = NA, dam = NA) %>%
  rbind(Families_first_brood %>%
          select(Dad) %>%
          rename(animal = Dad) %>%
          unique() %>%
          mutate(sire = NA, dam = NA)) %>%
  rbind(Families_first_brood %>%
          select(Offspring, Mom, Dad) %>%
          rename(animal = Offspring,
                 sire = Dad,
                 dam = Mom)) %>% 
  # Here we remove the offspring produced by female E4961 and that was not sequenced
  filter(!is.na(animal)) %>% 
  mutate(across(everything(), ~ as.integer(str_remove_all(., "E")))) %>% 
  as.data.frame()


########################################
########### Run animal model ###########
########################################
# traits <- c("Size", "Sum_curved_setaeP1P5_trunc", "Sum_spinesP6P7")
# 
# models_recap <- lapply(c("forsmani", "praehirsuta"), function(species, metadata, pedigree, traits){
#   print(species)
#   lapply(traits, function(trait, sp, metadata, pedigree){
#     print(paste0("    - Trait: ", trait))
# 
#     if (trait == "Size"){
#       effect_traits <- data.frame(
#         "Fixed_effects" = c(NA, "Sex"),
#         "Random_effects" = c(NA, NA),
#         "Model_ID" = c(1, 2)
#       )
#     }else if (trait == "Sum_curved_setaeP1P5_trunc"){
#       effect_traits <- data.frame(
#         "Fixed_effects" = c(NA, "Size", "Sum_spinesP6P7", "Size", "Sum_spinesP6P7"),
#         "Random_effects" = c(NA, NA, NA, NA, NA),
#         "Model_ID" = c(1, 2, 3, 4, 4)
#       )
#     } else {
#       effect_traits <- data.frame(
#         "Fixed_effects" = c(NA, "Size", "Sum_curved_setaeP1P5_trunc", "Size", "Sum_curved_setaeP1P5_trunc"),
#         "Random_effects" = c(NA, NA, NA, NA, NA),
#         "Model_ID" = c(1, 2, 3, 4, 4)
#       )
# 
#     }
# 
#     nb_cores_to_use <- min(detectCores() - 1, nrow(effect_traits))
#     cl <- makeCluster(nb_cores_to_use)
#     clusterEvalQ(cl, expr = {
#       library(tidyverse)
#       library(MCMCglmm)
#       })
# 
#     df_out <- parLapply(cl = cl,
#                         X = unique(effect_traits$Model_ID),
#                         fun = run_model_on_modelid,
#                         effect_traits, trait, sp, metadata, pedigree) %>%
#       bind_rows()
# 
#     stopCluster(cl)
#     return(df_out)
#   },
#   species, metadata, pedigree) %>%
#     bind_rows() %>%
#     return()
# }, metadata, pedigree, traits) %>%
#   bind_rows()
# 
# save(models_recap, file = "Associated_data/Rdata/Animal_models_recap.rdata")
load("Associated_data/Rdata/Animal_models_recap.rdata")
models_recap

########################################
########### Select best model ##########
########################################
# The model selection is done on the best DIC value (lowest)
# Here, the model choice is done on a combination of parameters: we try to select
# the model with the lowest DIC, if the autocorrelation is low. To look at this,
# we look at the effective sample sizes for the fixed and random parameters and
# select the model with the lowest DIC and the highest overall effective sample 
# sizes. For example, for the size of the individuals, for forsmani, both models
# converged (effective sample sizes above 900 for all parameters), so we chose the
# model with the lowest DIC. However, for the number of curved setae in forsmani,
# the model that converged best is the model number 4, but it is not the one with
# the lowest DIC (although this difference of DIC is not significant).

best_models <- data.frame(
  "Model_ID" = c(2, 2, 4, 4,
                 2, 3),
  "Trait" = c("Size", "Size", "Sum_curved_setaeP1P5_trunc", "Sum_curved_setaeP1P5_trunc",
              "Sum_spinesP6P7", "Sum_spinesP6P7"),
  "Species" = c("forsmani", "praehirsuta", "forsmani", "praehirsuta",
                "forsmani", "praehirsuta")
) %>% 
  left_join(models_recap, by = c("Model_ID", "Trait", "Species"))

########################################
########### Run best models ############
########################################
# best_models_recap <- lapply(unique(best_models$Species), function(species, metadata, pedigree, best_models){
#   print(species)
#   lapply(unique(best_models$Trait), function(trait, sp, metadata, pedigree, best_models){
#     print(paste0("    - Trait: ", trait))
#     #TODO: reprendre ça pour ne prendre que ce qu'il y a dans le tableau récap
#     effect_traits <- best_models %>%
#       filter(Trait == trait, Species == species) %>%
#       mutate(Fixed_effects = str_split(Effective_size_fixed_effects, "; ") %>%
#                lapply(function(x){str_remove_all(x, " [0-9\\.]+$")}) %>%
#                lapply(function(x){x[-c(x == "Intercept")]}) %>%
#                lapply(function(x){str_replace_all(x, "SexM", "Sex")}) %>%
#                lapply(function(x){if (length(x) == 0) return(NA) else return(x)}),
#              Random_effects = str_split(Effective_size_random_effects, "; ") %>%
#                lapply(function(x){str_remove_all(x, " [0-9\\.]+$")}) %>%
#                lapply(function(x){x[-c(x == "pedigree")]}) %>%
#                lapply(function(x){if (length(x) == 0) return(NA) else return(x)}))
# 
#     df_out <- lapply(unique(effect_traits$Model_ID), run_model_on_modelid,
#                         effect_traits, trait, sp, metadata, pedigree, TRUE) %>%
#       bind_rows()
# 
#     return(df_out)
#   },
#   species, metadata, pedigree, best_models) %>%
#     bind_rows() %>%
#     return()
# }, metadata, pedigree, best_models) %>%
#   bind_rows()
# 
# save(best_models_recap, file = "Associated_data/Rdata/Best_Animal_models_recap.rdata")
load("Associated_data/Rdata/Best_Animal_models_recap.rdata")

###################################################################
######### Look at posterior distributions of heritability #########
###################################################################
model <- readRDS("Associated_data/Heritability/data/Model_Sum_spinesP6P7_praehirsuta_3.rds")
summary(model)

heidel.diag(model[["Sol"]])


# Heritability
herit <- model$VCV[, "animal"] / rowSums(model$VCV)
effectiveSize(herit)
mean(herit)
HPDinterval(herit)
# plot(herit)

plot_list <- vector(mode = "list")
plot_list[[length(plot_list) + 1]] <- herit %>%
  as_tibble() %>% 
  rename(heritability = var1) %>% 
  rownames_to_column("Iteration") %>%
  mutate(Iteration = as.numeric(Iteration)) %>% 
  ggplot(aes(x = Iteration, y = heritability)) +
  geom_line() +
  my_theme

plot_list[[length(plot_list) + 1]] <- herit %>%
  as_tibble() %>% 
  rename(heritability = var1) %>% 
  ggplot(aes(x = heritability)) +
  geom_density(lwd = 1) +
  my_theme
for (colName in colnames(model$VCV)){
  plot_list[[length(plot_list) + 1]] <- model$VCV %>%
    as_tibble() %>%
    rownames_to_column("Iteration") %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    ggplot(aes(x = Iteration, y = !!sym(colName))) +
    geom_line() +
    my_theme
  
  plot_list[[length(plot_list) + 1]] <- model$VCV %>%
    as_tibble() %>%
    rownames_to_column("Iteration") %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    ggplot(aes(x = !!sym(colName))) +
    geom_density(lwd = 1) +
    my_theme
}

for (colName in colnames(model$Sol)){
  if (colName == "(Intercept)"){
    colName <- "Intercept"
  }
  plot_list[[length(plot_list) + 1]] <- model$Sol %>%
    as_tibble() %>%
    rownames_to_column("Iteration") %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    rename(Intercept = "(Intercept)") %>% 
    ggplot(aes(x = Iteration, y = !!sym(colName))) +
    geom_line() +
    my_theme
  
  plot_list[[length(plot_list) + 1]] <- model$Sol %>%
    as_tibble() %>%
    rownames_to_column("Iteration") %>%
    mutate(Iteration = as.numeric(Iteration)) %>%
    rename(Intercept = "(Intercept)") %>% 
    ggplot(aes(x = !!sym(colName))) +
    geom_density(lwd = 1) +
    my_theme
}






figure <- wrap_plots(plot_list) +
  plot_annotation(tag_level = "A") +
  plot_layout(design = "
  ACD
  BEF
  BGH
              ",
              widths = c(0.7, 0.3, 0.3),
              tag_level = "new")

figure
ggsave(plot = figure,
       "../Figures/Suppfigs/Posterior_distrib_size_praehirsuta_males.png",
       units = "px", scale = 4, height = 1200, width = 1400)


