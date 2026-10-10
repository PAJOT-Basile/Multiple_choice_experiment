#' Script to merge the metadata of the parents with that of the offspring (make the
#' tables compatible with one another).

#' This script is separated in two parts. First, merge the raw metadata, and then
#' add the species and verify the sex of the individuals (especially for some
#' immature individuals).

#' The first part of the script only uses the raw tsv/excel files that are provided
#' in the "Associated_data/Raw_data" folder.

#' The second part of the script requires data produced by genetic analysis (PCA
#' and parentage analysis), which require the first part of the script to be run
#' first. 

#' A good order to run this script is thus to run it until you save a file called
#' "raw_metadata.tsv", run scripts in the subfolder "Genomic_analysis" in this folder,
#' and then rerun this script until the end.


##############################
########## Libraries #########
##############################
libraries <- c("tidyverse", "readxl")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
########### Useful ###########
##############################
source("98_Transversal_scripts/Transversal_functions.r")

columns_in_first_df_not_in_second <- function(df1, df2){
  names_not_in_common <- names(df1)[which(names(df1) %not_in% names(df2))]
  return(names_not_in_common)
}
##############################
######## Import data #########
##############################
# Table of all individuals that were sequenced
# This table is important to add the individual id that is present in the VCF file
table_sequencing <- read_excel("Associated_data/Raw_metadata/Table_ddRAD_Manip_SexualSelection1_2025.xlsx",
                               sheet = 4,
                               col_names = TRUE) %>% 
  # We only keep the id of the individuals that is in the metadata tables and 
  # the id of the individuals that are in the VCF file
  select(ID_DNA_RAD, ID_ind0iv)



# Table of female parent individuals
table_females <- read.table("Associated_data/Raw_metadata/Female_parents.tsv",
                            sep = "\t", header = TRUE, dec = ",") %>% 
  # Here, we add a column that is invariant in this table, but will help us differentiate
  # parent individuals (that were put in the aquarium) from offspring.
  mutate(Family_level = "parents")



# Table of male parent individuals
table_males <- read.table("Associated_data/Raw_metadata/Father_phenotypes.tsv",
                            sep = "\t", header = TRUE, dec = ",") %>% 
  # In the same way as above, for females, we add an invariant column that will
  # later help us differentiate parents from offspring
  mutate(Family_level = "parents")



# Table of offspring phenotypes
table_offspring <- read.table("Associated_data/Raw_metadata/Offspring_phenotypes.tsv",
                              sep = "\t", header = TRUE, dec = ",") %>% 
  # In the same was as above, for the parents, we add an invariant column that
  # will help us differentiate offspring from the parents. Here the value is 
  # different from the ones above as this table contains the offspring data
  mutate(Family_level = "offspring",
  # We replace the name "hybrid" phenotype that can be misleading by the 
  # name "intermediate" phenotype, which is more appropriate
  Phenotype = ifelse(Phenotype == "hybrid", "Intermediate", Phenotype))


##############################
######## Merge data ##########
##############################
# The offspring table has the most columns, so we will start by reducing the number
# of columns in this table to match the columns in the other tables
table_offspring <- table_offspring %>%
  # First, for male parents, there is no spine classified as "very.large", 
  # while it is the case in offspring. So, we combine "large" and "very.large"
  # spines for all pereiopods (P7, P6, P5). The values of the "very.large" spines
  # will be stored in the columns "very.large.offspring" corresponding
  ## First, we rename the columns
  rename_with(~ paste0(., ".offspring"), c(P7.very.large.ep, P6.very.large.ep)) %>% 
  ## Then, we combine them
  mutate(P7.larg.ep = P7.larg.ep + P7.very.large.ep.offspring,
         P6.larg.ep = P6.larg.ep + P6.very.large.ep.offspring,
  # Then, for male parents, the detail of the size of the distant spines are not
  # given. Only the number of distant spines is given for P7. To merge them, we
  # compute the sum of the distant spines on P7. For the other columns, the
  # values are removed
  P7.ep.dist = P7.ep.dist.very.large + P7.ep.dist.large + P7.ep.dist.medium + P7.ep.dist.small,
  # We also recompute the sum of the number of spines on each pereiopd (P7-P4)
  P7.ep = P7.larg.ep + P7.med.ep + P7.small.ep,
  P6.ep = P6.larg.ep + P6.med.ep + P6.small.ep) %>% 
  # Once this is done, we remove the columns containing the detail of the size of distant spines
  select(-c(P7.ep.dist.very.large, P7.ep.dist.large, P7.ep.dist.medium, P7.ep.dist.small,
            P6.ep.dist.very.large, P6.ep.dist.large, P6.ep.dist.medium, P6.ep.dist.small,
            P5.ep.dist)) %>% 
  # We also remove all the columns that contain information on the ischium of individuals
  # and the number of spines on the P4 an P5 as there is a very small variability on these 
  # traits and they were not looked at for the parents
  select(-c(contains("isch"), P4.small.ep, P5.very.large.ep, P5.large.ep,
            P5.medium.ep, P5.small.ep, P5.ep)) %>% 
  # We also remove columns that contain the total number of curved setae on a pereiopod
  # as they need to be re-calculated with parents and offspring together (due to
  # the ">" that we need to study with parents and offspring together).
  select(-c(P5.curv.setae, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae)) %>% 
  # We also remove columns that were useful to follow the library preparation/sequencing 
  # process while it was going on
  select(-c(Sequencing, Library_prepared, Extraction_plate))





# Now, we need to prepare the table of the male parents to make it compatible with
# the table of the offspring
table_males <- table_males %>%
  # First, the names of the columns in the parent's table are not written in the 
  # same way as the names of the columns in the offspring's table, so we correct
  # this (for the columns Label, Phenotype, Size and Alcool_ID)
  rename_with(.cols = c(label, phenotype, size), ~ str_to_title(.x)) %>% 
  rename(Alcool_ID = alcool_ID) %>% 
  # Then, we add the columns that are missing in this table:
  #   The id of the mother: we do not know who the mother of the parents are when taken out of the aquarium
  mutate(Mother_ID = NA,
  #   Nor do we know when the parents were isolated
         Isolation_date = NA,
  # All the male parents are males
         Sex = "M",
  # As we transformed this column in the offspring's table, the number of "very.large" spines
  # for P5-P7 are not counted for the adults, so we add the column "very.large.offspring"
  # for the pereiopods and count them as NA (they were not counted)
  P7.very.large.ep.offspring = NA,
  P6.very.large.ep.offspring = NA) %>% 
  # Then, we rename the column containing the total number of distant spines on P7
  rename(P7.ep.dist = P7.ep.mid) %>% 
  # We add the sum of the number of spines on P7 and P6
  mutate(P7.ep = P7.larg.ep + P7.med.ep + P7.small.ep,
         P6.ep = P6.larg.ep + P6.med.ep + P6.small.ep,
  # The detail of the number of curved setae on the carp of P6, so we add it with NA values
  P6.carp.curv.setae = NA,
  # We also add a "Note" column to keep that of the offspring table
  Note = NA) %>% 
  # Finally, we remove the columns containing the total number of curved setae
  # per pereiopod, just as it was done in the offspring's table as we will need to
  # recompute them after combining parents and offspring
  select(-c(P5.curv.setae, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae))






# Now, we need to prepare the table of the female parents to be able to bind it
# with the two other tables (offspring and male parents)
table_females <- table_females %>%
  # First, we rename the columns that are in the female parent's table but written
  # differently from the corresponding columns in the offspring table
  rename(Alcool_ID = alcool_ID,
         Label = label,
         Death_date = death,
         Size = size,
         Photo_ID = Photo_IDs,
         Note = Notes) %>% 
  # Then, we add missing columns. We can start with the id of the mother and the 
  # date at which the female parents were isolated. These are NA as we do not 
  # have this information for individuals that were taken out of the aquarium
  mutate(Mother_ID = NA,
         Isolation_date = NA,
  # These individuals are all females
  Sex = "F",
  # For females, contrary to males, we cannot distinguish the species on the phenotype,
  # so this column is NA as well
  Phenotype = NA,
  # For female parents, the detailed phenotype (number of spines for females) was
  # not given, so for the spines, we add all the columns that are present in the
  # other tables, but they have a value of NA as they were not counted
  ## First, for P7
  P7.very.large.ep.offspring = NA,
  P7.larg.ep = NA,
  P7.med.ep = NA,
  P7.small.ep = NA,
  P7.ep = P7.larg.ep + P7.med.ep + P7.small.ep,
  P7.ep.dist = NA,
  ## Then for P6
  P6.very.large.ep.offspring = NA,
  P6.larg.ep = NA,
  P6.med.ep = NA,
  P6.small.ep = NA,
  P6.ep = P6.larg.ep + P6.med.ep + P6.small.ep,
  ## We also add the number of curved setae on P6, which is 0 for females
  P6.carp.curv.setae = 0,
  ## Same for P5, we add the values for spines (NA) and curved setae (0)
  P5.prop.curv.setae = 0,
  P5.carp.curv.setae = 0,
  P5.merus.curv.setae = 0,
  ## P4
  P4.prop.curv.setae = 0,
  P4.carp.curv.setae = 0,
  P4.merus.curv.setae = 0,
  ## P3
  P3.prop.curv.setae = 0,
  P3.carp.curv.setae = 0,
  P3.merus.curv.setae = 0,
  ## P2
  P2.prop.curv.setae = 0,
  P2.carp.curv.setae = 0,
  P2.merus.curv.setae = 0,
  ## P1
  P1.prop.curv.setae = 0,
  P1.carp.curv.setae = 0,
  P1.merus.curv.setae = 0)





# Now, there are some columns (see underneath) that are in the females table that
# are not in the other tables. These columns will be added when we combine the tables
# as to have less repetitive code
##### Columns that are not in the male parents and offspring tables, but that are in the female parent's table
columns_in_first_df_not_in_second(table_females, table_males)
columns_in_first_df_not_in_second(table_females, table_offspring)



# The next step is to combine the tables together
full_table <- table_offspring %>% 
  # We start by adding the male's table to the offspring's table as both the
  # offspring and the male's table are missing columns in common in comparison 
  # to the female's table
  rbind(table_males) %>% 
  # We add information that is given in the female's table on the broods that were
  # produced by the females (number of offspring in a brood and the names of the
  # corresponding plates in which the offspring were reared). At this point, the
  # parentage analysis was not run, so we do not have this information for the 
  # male parents.
  mutate(nb_tot_juv_brood1 = NA,
         nb_juv_brood1 = NA,
         nb_tot_juv_brood2 = NA,
         nb_juv_brood2 = NA,
         Name_offsprings_brood1 = NA,
         Name_offsprings_brood2 = NA) %>% 
  # Once we added the columns, we bind the female parent's table
  rbind(table_females) %>% 
  # Then, we convert the type of the column "Alcool_ID" to be able to add the 
  # sequencing id of the individuals that were sequenced
  mutate(Alcool_ID = as.character(Alcool_ID)) %>% 
  left_join(table_sequencing %>% 
               mutate(ID_ind0iv = str_split_fixed(ID_ind0iv, "r", 2)[, 1]),
            by = join_by("Alcool_ID" == "ID_ind0iv")) %>% 
  # Here, we also filter the table to remove individuals whose extraction ID contains
  # "_0" as these are duplicates that were not sequenced
  filter(!grepl("_0", ID_DNA_RAD))

# Then, we calculate the sum of the number of setae, for each pereiopod and the 
# sum per individual. To do that, we need to make sure that the values are comparable.
# The first step is to look at values that have a ">" and take the minimum for each
# segment of each pereiopod
minimum_distiguishable_counts <- full_table %>% 
  select(contains("curv")) %>% 
  mutate(across(everything(), ~ as.character(.))) %>% 
  pivot_longer(everything(), names_to = "Trait", values_to = "Count") %>% 
  filter(grepl(">", Count)) %>% 
  mutate(Count = as.numeric(str_remove_all(Count, ">"))) %>% 
  group_by(Trait) %>% 
  summarize(min_distinguisable_count = min(Count)) %>% 
  ungroup()

# Once this is done, we need to modify the raw data to truncate the data (when needed)
full_table <- full_table %>% 
  mutate(across(contains("curv"), ~ as.character(.))) %>% 
  pivot_longer(contains("curv"), names_to = "Trait", values_to = "Count") %>% 
  mutate(Count = as.numeric(str_remove_all(Count, ">"))) %>% 
  left_join(minimum_distiguishable_counts, by = "Trait") %>% 
  mutate(Count_truncated = case_when(!is.na(min_distinguisable_count) & Count > min_distinguisable_count ~ min_distinguisable_count,
                           TRUE ~ Count)) %>% 
  select(-min_distinguisable_count) %>% 
  pivot_wider(names_from = Trait, values_from = contains("Count", ignore.case = FALSE)) %>%
  rename_with( ~ str_replace_all(., "Count_truncated", "trunc"), .cols = starts_with("Count_truncated")) %>% 
  rename_with( ~ str_replace_all(., "Count_", "default_"), .cols = starts_with("Count_")) %>% 
  # Then, we calculate the sum of the setae in the data
  rowwise() %>% 
  # For P1
  mutate(default_P1.curv.setae = sum(default_P1.prop.curv.setae, default_P1.carp.curv.setae, default_P1.merus.curv.setae),
    trunc_P1.curv.setae = sum(trunc_P1.prop.curv.setae, trunc_P1.carp.curv.setae, trunc_P1.merus.curv.setae),
    # For P2
    default_P2.curv.setae = sum(default_P2.prop.curv.setae, default_P2.carp.curv.setae, default_P2.merus.curv.setae),
    trunc_P2.curv.setae = sum(trunc_P2.prop.curv.setae, trunc_P2.carp.curv.setae, trunc_P2.merus.curv.setae),
    # For P3
    default_P3.curv.setae = sum(default_P3.prop.curv.setae, default_P3.carp.curv.setae, default_P3.merus.curv.setae),
    trunc_P3.curv.setae = sum(trunc_P3.prop.curv.setae, trunc_P3.carp.curv.setae, trunc_P3.merus.curv.setae),
    # For P4
    default_P4.curv.setae = sum(default_P4.prop.curv.setae, default_P4.carp.curv.setae, default_P4.merus.curv.setae),
    trunc_P4.curv.setae = sum(trunc_P4.prop.curv.setae, trunc_P4.carp.curv.setae, trunc_P4.merus.curv.setae),
    # For P5
    default_P5.curv.setae = sum(default_P5.prop.curv.setae, default_P5.carp.curv.setae, default_P5.merus.curv.setae),
    trunc_P5.curv.setae = sum(trunc_P5.prop.curv.setae, trunc_P5.carp.curv.setae, trunc_P5.merus.curv.setae),
    # For P6
    default_P6.curv.setae = sum(default_P6.carp.curv.setae),
    trunc_P6.curv.setae = sum(trunc_P6.carp.curv.setae),
    # Get the total number of spines and setae
    Sum_curved_setaeP1P5_default = sum(default_P1.curv.setae, default_P2.curv.setae, default_P3.curv.setae,
                                       default_P4.curv.setae, default_P5.curv.setae),
    Sum_curved_setaeP1P5_trunc = sum(trunc_P1.curv.setae, trunc_P2.curv.setae, trunc_P3.curv.setae,
                                       trunc_P4.curv.setae, trunc_P5.curv.setae),
    Sum_spinesP6P7 = sum(P6.ep, P7.ep)
         )


# Write the table to be used later !!
write.table(full_table, "Associated_data/raw_metadata.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)






##############################
###### Use genomic data ######
##############################
# This is the second part of the script. In this part, make sure that you have the
# results of the species names using the genome-wide PCA as well as the verified
# sexes using the genomic data. 
# Using script Genomic_analysis/01_Genomic_cleaning_and_genomic_PCA_determine_species.r,
# we found there is a very good correspondence between the phenotype of the male
# individuals, the genomic clustering and the species of the offspring defined as
# a combination of the parent's species. We thus choose to use the species of the
# parents to identify the species of the offspring, which is what is done underneath.


# Import the correspondence between the names of the individuals and their species
# (genome-wide PCA output) for the parents only.
load("Associated_data/Rdata/Correspondence_name_species.rdata")

# For the offspring, the name of the species requires the output of the parentage
# analysis using the species of both the parents. To do this, we import the output 
# of COLONY.
Families <- read.table("Associated_data/Colony_output/small_data.ParentPair",
                       sep = ",", header = TRUE,
                       comment.char = "!") %>% 
  select(-Probability) %>% 
  # Here, we add the information for a mother that was not genotyped (as we have
  # the name of the mother for offspring)
  mutate(InferredMum = ifelse(InferredMum == "#", "E4933", InferredMum),
  # And we correct the names of the mothers for offspring that were changed from
  # one plate to another during their rearing and for which COLONY identified
  # a discrepancy between the mother of a brood and the mother of each individual.
         InferredMum = case_when(OffspringID %in% c("E5447", "E5449") ~ "E4936",
                                 OffspringID == "E5453" ~ "E4935",
                                 TRUE ~ InferredMum)) %>% 
  # Here, we change the name of the individuals to be compatible with the names
  # that were given to the individuals in the PCA output
  mutate(across(everything(), ~ str_split_fixed(., "r|_|b|a", 2)[, 1])) %>% 
  unique() %>% 
  # We add the species of the father
  left_join(Species_names_parents,
            by = join_by("InferredDad" == "ID_DNA_RAD")) %>% 
  # We add the species of the mother
  left_join(Species_names_parents,
            by = join_by("InferredMum" == "ID_DNA_RAD"),
            suffix = c("_dad", "_mom")) %>% 
  # We add the species of two females that were not phenotyped, but produced offspring
  # (this was done by counting the individuals of each sex and species that were
  # put in the aquarium and counting the individuals of each sex and species
  # that were taken out of the aquarium).
  mutate(Species_mom = case_when(InferredMum == "E4933" ~ "forsmani",
                                 InferredMum == "E4958" ~ "praehirsuta",
                                 TRUE ~ Species_mom),
  # And we determine the species of the offspring using the species of the parents
         Species = case_when(Species_mom == Species_dad ~ Species_mom,
                             TRUE ~ "Hybrid")) %>% 
  rename(ID_DNA_RAD = OffspringID) %>% 
  select(ID_DNA_RAD, Species)




full_table <- full_table %>% 
  mutate(id_dna_rad = str_split_fixed(ID_DNA_RAD, "r|b|_|a", 2)[, 1]) %>% 
  # Add the information of the species of the individuals using the PCA output and
  # the output of the parentage analysis (for the offspring)
  left_join(Species_names_parents %>% 
              rbind(Families),
            by = join_by("id_dna_rad" == "ID_DNA_RAD")) %>% 
  # Add species for individuals that were not in PCA (parents)
  mutate(Species = case_when(Alcool_ID == "4706" ~ Phenotype,
                             Alcool_ID == "4958" ~ "praehirsuta",
                             is.na(Species) & Sex == "M" & Alcool_ID != "5277" ~ Phenotype,
  # This individual had an intermediate phenotype, but it was fixed dead and not
  # genotyped, so this individual's species cannot be determined by its phenotype.
                             Alcool_ID == "5277" ~ NA,
                             TRUE ~ Species),
  # Finally, we convert the size column to numeric
  Size = as.numeric(Size)) %>% 
  select(-id_dna_rad)

# Write the table
write.table(full_table, "Associated_data/metadata.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

