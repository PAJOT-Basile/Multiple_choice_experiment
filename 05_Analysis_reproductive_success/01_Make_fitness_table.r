#' This script uses the output of the parentage analysis to make one large table
#' that summarizes the reproductive success of each parent into several fitness
#' components that will be used in the later analysis.

###########################################################
######################## Libraries ########################
###########################################################
libraries <- c("tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)

###########################################################
##################### Useful functions ####################
###########################################################
source("98_Transversal_scripts/Transversal_functions.r")

###########################################################
####################### Import data #######################
###########################################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b")) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD))

# Import the outputs of colony
Families <- read.table("Associated_data/Colony_output/small_data.ParentPair",
                       sep = ",", header = TRUE, comment.char = "!") %>% 
  select(-Probability) %>% 
  # Remove the replicates
  mutate(across(everything(), ~ str_split_fixed(., "r", 2)[, 1] %>% 
                  str_remove_all("a|b"))) %>% 
  unique() %>% 
  # Here, we add a mom that was not genotyped as it was missing, but for which we 
  # know the offspring.
  mutate(InferredMum = case_when(InferredMum == "#" ~ "E4933",
                                 TRUE ~ InferredMum)) %>% 
  # Rename the columns so they are more user-friendly :)
  rename_with(~ c("Offspring", "Dad", "Mom"), .cols = everything())


# Import the survival per plate
survival <- read.table("Associated_data/Raw_metadata/Survival_offspring.tsv",
                       header = TRUE, sep = "\t") %>% 
  # Here, we need to have the survival of the whole brood, so to summarize the
  # number of offspring that are alive after n days in the whole brood
  group_by(Mother_ID, Isolation_date) %>% 
  summarize(Plate_ID = paste(Plate_ID, collapse = "-"),
            Nb_indivs_isolation = sum(Nb_indivs_isolation, na.rm = TRUE),
            Nb_indivs_40d = sum(Nb_indivs_40d, na.rm = TRUE),
            Nb_indivs_60d = sum(Nb_indivs_60d, na.rm = TRUE)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  # And to link this table with others, we need to have the name of the brood, which
  # in other tables is under the format: "Ofirst_plate-Olast_plate", so we get
  # get the names of the first and last brood and paste them to get the brood name
  mutate(First_plate = str_split_fixed(Plate_ID, "-", 4)[, 1],
         Last_plate = str_split_last(Plate_ID, "-"),
         Brood_name = ifelse(First_plate == Last_plate,
                             First_plate,
                             paste(First_plate, Last_plate, sep = "-"))) %>% 
  # And we compute the survival per brood using the number of offspring alive at
  # n days and divide it by the number of offspring that were isolated
  mutate(Surv_40d = Nb_indivs_40d / Nb_indivs_isolation,
         Surv_60d = Nb_indivs_60d / Nb_indivs_isolation) %>% 
  # Finally, we clean the table
  select(-c(Isolation_date, Plate_ID, Nb_indivs_isolation, Nb_indivs_40d, Nb_indivs_60d, First_plate, Last_plate, Mother_ID))


##########################################################
############# Make table of fitness components ###########
##########################################################
# For each parent, we would like to have several fitness components:
#   1. The number of partners they had
#   2. The number of offspring they parented
#   3. The contribution of each parent to the produced broods
#   4. The survival of all the offspring of each parent

# To do this, we need to add to the Family table that we imported the number of 
# offspring that were produced in the first brood as well as the survival of the
# each brood.
Families_first_brood <- Families %>%
  # Here, we add the name of the brood for each offspring to only keep the offspring
  # that were produced in the first brood.
  left_join(metadata %>% 
              select(ID_DNA_RAD, Name_offsprings_brood1, nb_tot_juv_brood1),
            by = join_by("Mom" == "ID_DNA_RAD")) %>% 
  mutate(Brood_name = Name_offsprings_brood1) %>% 
  # As the name of the broods are under the format: "Oplate_nb1-Oplate_nb2",
  # we separate this name to keep only the number of the plates, so as to keep
  # individuals in these plates.
  separate(Name_offsprings_brood1, c("min_brood_nb", "max_brood_nb"), sep = "-") %>% 
  mutate(across(ends_with("brood_nb"), ~ str_remove_all(., "O") %>% 
                  as.numeric()),
  # In the case where only one plate was isolated, the second brood number is not
  # defined, so we replace this with the number of the plate
         max_brood_nb = ifelse(is.na(max_brood_nb), min_brood_nb, max_brood_nb)) %>% 
  # Then, we add the label of the offspring, that is under the following format:
  # "Oplate_nb-well_nb". We only want to keep the number of the plate to keep
  # individuals that are in the plate numbers identified earlier.
  left_join(metadata %>% 
              select(ID_DNA_RAD, Label) %>% 
              mutate(Plate_nb = str_split_fixed(Label, "-", 2)[, 1] %>% 
                       str_remove_all("O") %>% 
                       as.numeric()),
            join_by("Offspring" == "ID_DNA_RAD")) %>% 
  # We then filter offspring that are in the first brood that was parented.
  filter(Plate_nb <= max_brood_nb & Plate_nb >= min_brood_nb) %>% 
  # And we clean up the table
  rename(nb_tot_juvs = nb_tot_juv_brood1) %>% 
  select(-c(ends_with("brood_nb"), Label, Plate_nb))


save(Families_first_brood, file = "Associated_data/Rdata/Families_first_brood.rdata")

# The next step is to add the contribution of each parent to the produced broods.
# To do this, we need to have a summary table containing the number of offspring
# that were born per brood and the number of offspring that were genotyped from
# each brood
Brood_sizes <- Families_first_brood %>% 
  group_by(Mom, Brood_name) %>% 
  summarize(nb_offspring_per_brood = unique(nb_tot_juvs),
            seq_indivs = n_distinct(Offspring))


# Then, we compute the contribution of the male parents to each brood (for the
# female parents, the contribution to the brood is always 1 as there cannot be
# several mothers in a brood). To do this, we divide the number of offspring sired
# by a male per brood by the number of offspring that were sequenced per brood
Contribution_of_dad_to_broods <- Families_first_brood %>% 
  # Here, we calculate the number of offspring sired by male parents per brood
  group_by(Dad, Brood_name) %>% 
  summarize(nb_offspring_sired_per_brood = n_distinct(Offspring)) %>% 
  # Then, using the name of the brood, we add the number of genotyped individuals
  # in the brood to compute the contribution
  left_join(Brood_sizes,
            by = "Brood_name") %>% 
  mutate(Contribution = nb_offspring_sired_per_brood / seq_indivs) %>% 
  select(Dad, Brood_name, Contribution)

  
# Finally, we need to combine the variables (number of offspring, contribution
# to a brood, survival of the offspring) and add the number of partners of the 
# parents into a single table. To do this, males and females will be treated
# separately. For males, the total number of sired offspring if the focal male was
# the only one contributing to a brood, the total number of sired offspring and
# the total number of offspring alive at n days are first calculated per brood
# and then summed.
fitness_components_males <- survival %>%
  # Here, we add the number of offspring per brood as well as the parent's ids
  left_join(Families_first_brood %>% 
              select(Brood_name, nb_tot_juvs, Mom, Dad),
            by = "Brood_name") %>% 
  # Here, we add the contribution of the male parents to the different broods
  left_join(Contribution_of_dad_to_broods,
            by = c("Brood_name", "Dad")) %>% 
  unique() %>% 
  # Here, we remove broods that are second broods, but that were not removed from
  # the data
  filter(!is.na(nb_tot_juvs)) %>% 
  # And we compute the different stats:
  #   1. The total number of sired offspring with the contribution of the male
  #       to the brood
  mutate(nb_offspring_produced_by_dad = Contribution * nb_tot_juvs,
  #   2. The total number of offspring alive at n days
         nb_offspring_produced_by_dad_alive_40 = Contribution * nb_tot_juvs * Surv_40d,
         nb_offspring_produced_by_dad_alive_60 = Contribution * nb_tot_juvs * Surv_60d) %>% 
  group_by(Dad) %>% 
  # Finally, we compute the different metrics per indivdiual
  #   1. The number of partners
  summarize(nb_partners = n_distinct(Mom),
  #   2. The total number of sired offspring if the focal male was the only one
  #       contributing to a brood
            max_nb_offspring = sum(nb_tot_juvs, na.rm = TRUE),
    # 3. The total number of sired offspring with the actual contribution of the
    #     focal male to the brood
            nb_offspring = sum(nb_offspring_produced_by_dad, na.rm = TRUE),
  #   4. The total number of sired offspring alive at n days
            nb_offspring_surv40 = sum(nb_offspring_produced_by_dad_alive_40, na.rm = TRUE),
            nb_offspring_surv60 = sum(nb_offspring_produced_by_dad_alive_60, na.rm = TRUE)) %>% 
  # In the previous table, all the male parents that were added to the aquarium are
  # not found as some of them did not have any offspring. However, the fitness of
  # these males does exist and it is 0 for all the above fitness components, so we
  # add them to the table.
  right_join(metadata %>% 
               filter(Family_level == "parents",
                      Sex == "M") %>% 
               select(ID_DNA_RAD, Species, Sex),
             by = join_by("Dad" == "ID_DNA_RAD")) %>% 
  # Here, we added the males that did not reproduce, but their fitness values
  # in the table are equal to NA. They however should be equal to 0, so we modify
  # them.
  mutate(across(!c(Dad, Species, Sex), ~ case_when(is.na(.) ~ 0,
                                  TRUE ~ .))) %>% 
  rename(Parent_ID = Dad)

# Once this is done, we do the same for the females
fitness_components_females <- survival %>%
  # Here, we add the number of offspring per brood as well as the parent's ids
  left_join(Families_first_brood %>% 
              select(Brood_name, nb_tot_juvs, Mom, Dad),
            by = "Brood_name") %>% 
  unique() %>% 
  # Here, we remove broods that are second broods, but that were not removed from
  # the data
  filter(!is.na(nb_tot_juvs)) %>% 
  # And we compute the different stats. We already have the total number of
  # offspring parented by each female, so the remaining stats are:
  #   The total number of sired offspring alive at n days
  mutate(nb_offspring_alive_40 = Surv_40d * nb_tot_juvs,
         nb_offspring_alive_60 = Surv_60d * nb_tot_juvs) %>% 
  group_by(Mom) %>% 
  # Finally, we compute the different metrics per indivdiual
  #   1. The number of partners
  summarize(nb_partners = n_distinct(Dad),
  #   2. The total number of sired offspring if the focal male was the only one
  #       contributing to a brood
            max_nb_offspring = unique(nb_tot_juvs),
  #   3. The total number of sired offspring with the actual contribution of the 
  #       focal male to the brood
  # (The two components calcualted above are the same for females)
            nb_offspring = sum(unique(nb_tot_juvs), na.rm = TRUE),
  #   4. The total number of sired offspring alive at n days
            nb_offspring_surv40 = sum(unique(nb_offspring_alive_40), na.rm = TRUE),
            nb_offspring_surv60 = sum(unique(nb_offspring_alive_60), na.rm = TRUE)) %>% 
  # In the previous table, all the male parents that were added to the aquarium are
  # not found as some of them did not have any offspring. However, the fitness of
  # these males does exist and it is 0 for all the above fitness components, so we
  # add them to the table.
  right_join(metadata %>% 
               filter(Family_level == "parents",
                      Sex == "F") %>% 
               select(ID_DNA_RAD, Species, Sex),
             by = join_by("Mom" == "ID_DNA_RAD")) %>% 
  # Here, we added the males that did not reproduce, but their fitness values
  # in the table are equal to NA. They however should be equal to 0, so we modify
  # them.
  mutate(across(!c(Mom, Species, Sex), ~ case_when(is.na(.) ~ 0,
                                  TRUE ~ .))) %>% 
  rename(Parent_ID = Mom) %>% 
  # Here, one female had an offspring that did not survive 40 days and was not
  # genotyped. We thus know this female reproduced with a male, but we do not
  # know which male it reproduced with. We thus chose to remove this female from
  # the analyses
  filter(Parent_ID != "E4961")
  


# Save the reproductive success of each individual
fitness_components <- fitness_components_males %>% 
  rbind(fitness_components_females) %>% 
  relocate(c(Species, Sex), .after = Parent_ID)

save(fitness_components, file = "Associated_data/Rdata/Fitness_components_parents.rdata")