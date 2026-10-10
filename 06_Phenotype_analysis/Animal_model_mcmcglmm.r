########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "MCMCglmm")
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
  mutate(across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()

# Make the pedigree
load("Associated_data/Rdata/Families_first_brood.rda")
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
  mutate(across(everything(), ~ as.integer(str_remove_all(., "E")))) %>% 
  as.data.frame()

metadata_sp <- metadata  %>%
  filter(animal %in% pedigree$animal,
         Species == "forsmani") %>% 
  as.data.frame()

########################################
########### Test simple model ##########
########################################
prior <- list(R = list(V = 1, nu = 0.002),
              G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

model <- MCMCglmm(Sum_spinesP6P7 ~ 1,
                  random = ~ animal,
                  family = "poisson",
                  prior = prior,
                  pedigree = pedigree,
                  data = metadata_sp,
                  nitt = 100000,
                  burnin = 10000,
                  thin = 10)

plot(model$Sol)
effectiveSize(model[["Sol"]])
autocorr.diag(model[["Sol"]])
heidel.diag(model[["Sol"]])


plot(model$VCV)
effectiveSize(model[["VCV"]])
autocorr.diag(model[["VCV"]])
heidel.diag(model[["VCV"]])

summary(model)


# Heritability
herit <- model$VCV[, "animal"] / rowSums(model$VCV)
effectiveSize(herit)
mean(herit)
HPDinterval(herit)
plot(herit)
