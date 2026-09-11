########################################
############## Libraries ###############
########################################
libraries <- c("tidyverse", "adegenet", "vcfR", "readxl", "statgenGWAS", "ggforce", "ggh4x", "ggpubr")
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
