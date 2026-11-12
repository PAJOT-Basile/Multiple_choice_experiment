########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "MCMCglmm", "parallel")
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
########### Test simple model ##########
########################################
traits <- c("Size", "Sum_curved_setaeP1P5_trunc", "Sum_spinesP6P7")

# TODO: if you have to rerun this part, try putting it in parallel to win some time

models_recap <- lapply(c("forsmani", "praehirsuta"), function(species, metadata, pedigree, traits){
  print(species)
  lapply(traits, function(trait, sp, metadata, pedigree){
    print(paste0("    - Trait: ", trait))

    if (trait == "Size"){
      effect_traits <- data.frame(
        "Fixed_effects" = c(NA, "Sex", NA),
        "Random_effects" = c(NA, NA, "Sex"),
        "Model_ID" = c(1, 2, 3)
      )
    }else if (trait == "Sum_curved_setaeP1P5_trunc"){
      effect_traits <- data.frame(
        "Fixed_effects" = c(NA, "Size", "Sum_spinesP6P7", "Size", "Sum_spinesP6P7",
                            NA, NA, NA, NA, "Size", "Sum_spinesP6P7"),
        "Random_effects" = c(NA, NA, NA, NA, NA, "Size", "Sum_spinesP6P7",
                             "Size", "Sum_spinesP6P7", "Sum_spinesP6P7", "Size"),
        "Model_ID" = c(1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9)
      )
    } else {
      effect_traits <- data.frame(
        "Fixed_effects" = c(NA, "Size", "Sum_curved_setaeP1P5_trunc", "Size", "Sum_curved_setaeP1P5_trunc",
                            NA, NA, NA, NA, "Size", "Sum_curved_setaeP1P5_trunc"),
        "Random_effects" = c(NA, NA, NA, NA, NA, "Size", "Sum_curved_setaeP1P5_trunc",
                             "Size", "Sum_curved_setaeP1P5_trunc", "Sum_curved_setaeP1P5_trunc", "Size"),
        "Model_ID" = c(1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 9)
      )

    }

    lapply(unique(effect_traits$Model_ID), function(model_id, effect_traits, trait, sp, metadata, pedigree){
      print(paste0("       - Model_ID: ", model_id))

      fixed_effect_traits <- effect_traits %>%
        filter(Model_ID == model_id) %>%
        pull(Fixed_effects)
      random_effect_traits <- effect_traits %>%
        filter(Model_ID == model_id) %>%
        pull(Random_effects)
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
               Species == species)

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
      variances_traits <- lapply(as.vector(na.omit(c(trait, random_effect_traits))), function(x, df){
        if (x == "Sex") return(NULL)
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
      priors_random <- lapply(as.vector(na.omit(c(trait, random_effect_traits))), function(x, df){
        if (x == "Sex") return(NULL)
        v <- df %>%
          filter(Trait == x) %>%
          pull(Variance)

        if (v >= 1){
          prior_val <- list(V = 1, nu = 0.002)
        } else {
          prior_val <- list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)
        }
        return(prior_val)
      }, variances_traits)

      names(priors_random) <- paste0("G", 1:length(as.vector(na.omit(c(trait, random_effect_traits)))))


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
                        # nitt = 2000000,
                        # burnin = 200000,
                        # thin = 2000)
                        nitt = 100,
                        burnin = 10,
                        thin = 10)


      # Get some output stats
      DIC_mod <- summary(model)$DIC
      effective_size_fixed_effects <- effectiveSize(model[["Sol"]])
      effective_size_random_effects <- effectiveSize(model[["VCV"]])[names(effectiveSize(model[["VCV"]])) != "units"]

      # Store the output in a data frame
      out_df <- data.frame("Model_ID" = model_id,
                           "Trait" = trait,
                           "Species" = species,
                           "DIC" = DIC_mod,
                           "Effective_size_fixed_effects" = paste(fixed_effect_traits, effective_size_fixed_effects,
                                                                  collapse = "; "),
                           "Effective_size_random_effects" = paste(random_effect_traits, effective_size_random_effects,
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

      return(out_df)
    }, effect_traits, trait, sp, metadata, pedigree) %>%
      bind_rows() %>%
      return()

  },
  species, metadata, pedigree) %>%
    bind_rows() %>%
    return()
}, metadata, pedigree, traits) %>%
  bind_rows()

# save(models_recap, file = "Associated_data/Rdata/Animal_models_recap.rdata")
load("Associated_data/Rdata/Animal_models_recap.rdata")











plot(model$Sol)
effectiveSize(model[["Sol"]])
# autocorr.diag(model[["Sol"]])
# heidel.diag(model[["Sol"]])


plot(model$VCV)
effectiveSize(model[["VCV"]])
# autocorr.diag(model[["VCV"]])
# heidel.diag(model[["VCV"]])

summary(model)


# Heritability
herit <- model$VCV[, "animal"] / rowSums(model$VCV)
effectiveSize(herit)
mean(herit)
HPDinterval(herit)
plot(herit)
