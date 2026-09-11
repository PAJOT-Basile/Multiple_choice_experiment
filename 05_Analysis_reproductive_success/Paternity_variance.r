########## Libraries ##########
libraries <- c("tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

########## Functions ##########
"%not_in%" <- function(x, y) return(!(x %in% y))

str_split_last <- function(string, pattern){
  nb_substrings <- str_count(string, pattern) + 1
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

########## Import data ##########
# Import the metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1],
         ID_DNA_RAD = str_remove_all(ID_DNA_RAD, paste(c("a", "b"), collapse = "|")),
         Tot_juvs = nb_tot_juv_brood1 + nb_tot_juv_brood2,
         Size = str_replace_all(Size, ",", ".") %>% 
           as.numeric()) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD))

# Import the outputs of colony
paternity_results <- read.table("/shared/projects/sexisol/Softwares/Colony/small_data_with_mothers/small_data.Paternity",
                                sep = ",", header = TRUE) %>% 
  # remove sequencing duplicates
  mutate(OffspringID = str_split_fixed(OffspringID, "r", 2)[, 1],
         OffspringID = str_remove_all(OffspringID, paste(c("a", "b"), collapse = "|"))) %>% 
  unique()




########## Prepare data ##########
# Make a table of offspring and mothers
known_maternity <- metadata %>% 
  filter(Family_level == "offspring") %>% 
  select(ID_DNA_RAD, Mother_ID) %>% 
  rename(Offspring = ID_DNA_RAD) %>% 
  left_join(metadata %>% 
              filter(Sex == "F", Family_level == "parents") %>% 
              select(Label, ID_DNA_RAD) %>% 
              rename(Mom = ID_DNA_RAD)  %>% 
              # Remove controls that
              filter(!grepl("_0", Mom)) %>% 
              # remove sequencing duplicates
              mutate(Mom = str_split_fixed(Mom, "r", 2)[, 1],
                     Mom = str_remove_all(Mom, paste(c("a", "b"), collapse = "|"))),
            by = join_by("Mother_ID" == "Label"),
            relationship = "many-to-one") %>% 
  select(-Mother_ID)

# Make a family table with the id of the offspring and both its parents
Families <- known_maternity %>% 
  # Add the dad
  left_join(paternity_results %>% 
              select(-ProbDad1) %>% 
              rename(Dad = InferredDad1),
            by = join_by("Offspring" == "OffspringID")) %>% 
  # Add the species of the parents
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species),
            by = join_by("Mom" == "ID_DNA_RAD")) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species),
            by = join_by("Dad" == "ID_DNA_RAD"),
            suffix = c("_mom", "_dad")) %>% 
  # Add species of offspring depending on the parents' species
  mutate(Species_offspring = case_when(
    Species_mom == "praehirsuta" & Species_dad == "praehirsuta" ~ "praehirsuta",
    Species_mom == "forsmani" & Species_dad == "forsmani" ~ "forsmani",
    is.na(Species_mom) | is.na(Species_dad) ~ NA,
    TRUE ~ paste0(
      str_split_fixed(Species_mom, "", 2)[, 1],
      str_split_fixed(Species_dad, "", 2)[, 1],
      "_hybrid"
    )
  )) %>% 
  group_by(Mom, Dad) %>%
  mutate(Family_ID = cur_group_id()) %>% 
  ungroup()

save(Families, file = "/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Sexual_selection/Families.rda")
########## Hybrids ##########
hybrid_offspring <- Families %>%
  filter(grepl("hybrid", Species_offspring))

print(paste("There are", nrow(hybrid_offspring), "hybrids in the dataset."))

nb_hybrids <- hybrid_offspring %>% 
  nrow()

########## Get the brood name ##########
Families <- Families %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Name_offsprings_brood1, Name_offsprings_brood2, nb_tot_juv_brood1, nb_tot_juv_brood2),
            by = join_by("Mom" == "ID_DNA_RAD")) %>% 
  mutate(min_brood_nb = str_split_fixed(Name_offsprings_brood1, "-", 2)[, 1] %>% 
      str_remove_all("O") %>% 
      as.numeric(),
    max_brood_nb = str_split_fixed(Name_offsprings_brood1, "-", 2)[, 2] %>% 
        str_remove_all("O") %>% 
        as.numeric(),
    max_brood_nb = ifelse(is.na(max_brood_nb), min_brood_nb, max_brood_nb)) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Label),
            join_by("Offspring" == "ID_DNA_RAD")) %>% 
  mutate(Plate = str_split_fixed(Label, "-", 2)[, 1] %>% 
           str_remove_all("O") %>% 
           as.numeric()) %>% 
  select(-c(Name_offsprings_brood2, nb_tot_juv_brood2)) %>% 
  rename(Brood_name = Name_offsprings_brood1,
         nb_tot_juvs = nb_tot_juv_brood1)

# Keep only individuals from the first brood
Families_first_brood_only <- Families %>% 
  # We filter individuals that are in the first brood and keep the two individuals indicated
  # at the end as they were individuals from the first brood that were not isolated in the right 
  # plates.
  filter(Plate <= max_brood_nb & Plate >= min_brood_nb | Offspring %in% c("E5111", "E5453")) %>% 
  select(-c(min_brood_nb, max_brood_nb, Label, Plate))

save(Families_first_brood_only, file = "/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Sexual_selection/Families.rda")


Brood_sizes <- Families_first_brood_only %>% 
  group_by(Mom, Brood_name) %>% 
  summarize(nb_offspring_per_brood = unique(nb_tot_juvs),
            seq_indivs = n_distinct(Offspring))


Contribution_of_dad_to_broods <- Families_first_brood_only %>% 
  group_by(Dad, Brood_name) %>% 
  summarize(nb_offspring = n_distinct(Offspring)) %>% 
  left_join(Brood_sizes,
            by = "Brood_name") %>% 
  mutate(Contribution_to_brood = round(nb_offspring / seq_indivs, digits = 2)) %>% 
  select(Dad, Brood_name, Contribution_to_brood)

survival <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/Survival_offspring.tsv",
                       header = TRUE, sep = "\t") %>% 
  group_by(Mother_ID, Isolation_date) %>% 
  summarize(Plate_ID = paste(Plate_ID, collapse = "-"),
            Nb_indivs_isolation = sum(Nb_indivs_isolation, na.rm = TRUE),
            Nb_indivs_40d = sum(Nb_indivs_40d, na.rm = TRUE),
            Nb_indivs_60d = sum(Nb_indivs_60d, na.rm = TRUE)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(First_plate = str_split_fixed(Plate_ID, "-", 4)[, 1],
         Last_plate = str_split_last(Plate_ID, "-"),
         Brood_name = ifelse(First_plate == Last_plate, First_plate, paste(First_plate, Last_plate, sep = "-"))) %>% 
  mutate(Surv_40d = Nb_indivs_40d / Nb_indivs_isolation,
         Surv_60d = Nb_indivs_60d / Nb_indivs_isolation) %>% 
  select(-c(Isolation_date, Plate_ID, Nb_indivs_isolation, Nb_indivs_40d, Nb_indivs_60d, First_plate, Last_plate, Mother_ID))
  


survival_offspring_dad <- survival %>%
  left_join(Families_first_brood_only %>% 
              select(Brood_name, nb_tot_juvs, Mom, Dad),
            by = "Brood_name") %>% 
  left_join(Contribution_of_dad_to_broods,
            by = c("Brood_name", "Dad")) %>% 
  unique() %>% 
  filter(!is.na(nb_tot_juvs)) %>% 
  mutate(nb_offspring_produced_by_dad = Contribution_to_brood * nb_tot_juvs,
         nb_offspring_produced_by_dad_alive_40 = Contribution_to_brood * nb_tot_juvs * Surv_40d,
         nb_offspring_produced_by_dad_alive_60 = Contribution_to_brood * nb_tot_juvs * Surv_60d) %>% 
  group_by(Dad) %>% 
  summarize(nb_mates = n_distinct(Mom),
            max_nb_offspring = sum(nb_tot_juvs, na.rm = TRUE),
            nb_offspring = sum(nb_offspring_produced_by_dad, na.rm = TRUE),
            nb_offspring_surv40 = sum(nb_offspring_produced_by_dad_alive_40, na.rm = TRUE),
            nb_offspring_surv60 = sum(nb_offspring_produced_by_dad_alive_60, na.rm = TRUE))



########## Paternity results ##########
# Look at the number of parents that reproduced
## fathers
Contribution_of_dad_to_broods %>% 
  ungroup() %>% 
  pull(Dad) %>% 
  unique() %>% 
  length()

all_fathers <- metadata %>% 
  filter(Family_level == "parents",
         Sex == "M") %>% 
  pull(ID_DNA_RAD)

## mothers
Brood_sizes %>% 
  pull(Mom) %>% 
  unique() %>% 
  length()

all_mothers <- metadata %>% 
  filter(Family_level == "parents",
         Sex == "F") %>% 
  pull(ID_DNA_RAD)


# Look at the reproductive success of parents
# Fathers
fathers_partners_and_offspring <- survival_offspring_dad %>% 
  rbind(data.frame(
    "Dad" = all_fathers[which(all_fathers %not_in% survival_offspring_dad$Dad)],
    "nb_mates" = 0,
    "max_nb_offspring" = 0,
    "nb_offspring" = 0,
    "nb_offspring_surv40" = 0,
    "nb_offspring_surv60" = 0
  )) %>% 
  mutate(Parent = "Males") %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species),
            by = join_by("Dad" == "ID_DNA_RAD")) %>% 
  rename(Parent_ID = Dad)

# Mothers
survival_offspring_mom <- survival %>%
  left_join(Families_first_brood_only %>% 
              select(Brood_name, nb_tot_juvs, Mom, Dad),
            by = "Brood_name") %>% 
  unique() %>% 
  filter(!is.na(nb_tot_juvs)) %>% 
  mutate(nb_offspring_alive_40 = Surv_40d * nb_tot_juvs,
         nb_offspring_alive_60 = Surv_60d * nb_tot_juvs) %>% 
  group_by(Mom) %>% 
  summarize(nb_mates = n_distinct(Dad),
            max_nb_offspring = unique(nb_tot_juvs),
            nb_offspring = sum(unique(nb_tot_juvs), na.rm = TRUE),
            nb_offspring_surv40 = sum(unique(nb_offspring_alive_40), na.rm = TRUE),
            nb_offspring_surv60 = sum(unique(nb_offspring_alive_60), na.rm = TRUE))

mothers_partners_and_offspring <- survival_offspring_mom %>% 
  rbind(data.frame(
    "Mom" = all_mothers[which(all_mothers %not_in% survival_offspring_mom$Mom)],
    "nb_mates" = 0,
    "max_nb_offspring" = 0,
    "nb_offspring" = 0,
    "nb_offspring_surv40" = 0,
    "nb_offspring_surv60" = 0
  )) %>% 
  mutate(Parent = "Females") %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Species),
            by = join_by("Mom" == "ID_DNA_RAD")) %>%
  rename(Parent_ID = Mom)

# Save the reproductive success of each individual
reproductive_success_parents <- fathers_partners_and_offspring %>% 
  rbind(mothers_partners_and_offspring)

save(reproductive_success_parents, file = "/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Sexual_selection/reproductive_success_parents.Rdata")

# Look at the variance of reproductive success
reproductive_success_parents %>% 
  group_by(Parent, Species) %>% 
  summarize(mean_nb_mates = mean(nb_mates),
            var_nb_mates = var(nb_mates),
            mean_max_nb_offspring = mean(max_nb_offspring),
            var_max_nb_offspring = var(max_nb_offspring),
            mean_nb_offspring = mean(nb_offspring),
            var_nb_offspring = var(nb_offspring),
            mean_nb_offspring_surv40 = mean(nb_offspring_surv40),
            var_nb_offspring_surv40 = var(nb_offspring_surv40),
            mean_nb_offspring_surv60 = mean(nb_offspring_surv60),
            var_nb_offspring_surv60 = var(nb_offspring_surv60))

