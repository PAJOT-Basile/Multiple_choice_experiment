########## Libraries ##########
libraries <- c("tidyverse")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)
rm(libraries)



########## Import data ##########
# Import the metadata
metadata <- read.table("Associated_data/raw_metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that were not added in the sequencing
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
  mutate(Tot_juvs = nb_tot_juv_brood1 + nb_tot_juv_brood2)

# Import the outputs of colony. We can only check if colony worked well using the mothers.
# We know which mothers gave which offspring as we have direct development, but
# we can also check this using colony
maternity_results <- read.table("/shared/projects/sexisol/Softwares/Colony/small_data_without_mothers/small_data_without_mothers.Maternity",
                                sep = ",", header = TRUE, fill = TRUE)



########## Maternity results ##########
# Look if the mothers were well assigned
maternity_comparison <- metadata %>% 
  filter(Family_level == "offspring") %>% 
  select(ID_DNA_RAD, Mother_ID) %>% 
  # Join to the mothers to get their ID
  left_join(metadata %>% 
              select(ID_DNA_RAD, Label),
            by = join_by("Mother_ID" == "Label"),
            suffix = c("", "_mother")) %>% 
    rename(Known_mother = ID_DNA_RAD_mother) %>% 
  # Add the maternity outputs
  left_join(maternity_results,
            by = join_by("ID_DNA_RAD" == "OffspringID")) %>% 
  # Remove individuals that were not attributed to a mother as the mothers
  # were not added to the genotyping
  filter(!is.na(InferredMum1)) %>% 
  # Check if the known mother is the inferred mother (we only check the first mother)
  mutate(Is_mother_well_infered = (Known_mother == InferredMum1))

# We say that the mother is badly inferred if the name of the inferred mother is 
# different from the name of the known mother.
# Here, we filter to keep only the mothers that are badly inferred
mothers_bad_assignment <- maternity_comparison %>% 
  filter(!Is_mother_well_infered)

# And we calculate the percentage of mothers that were badly inferred on the 
# overall dataset.
perc_bad_infer <- round((nrow(mothers_bad_assignment) / nrow(maternity_comparison)) * 100,
                        digits = 2)

print(paste(perc_bad_infer, "% of maternal relationships were badly inferred."))

# Check to see if there is a correlation with the brood size (if smaller broods
# are less well attributed)
relationship_brood_size_bad_attribution <- maternity_comparison %>% 
  group_by(Known_mother) %>% 
  summarize(count = n()) %>% 
  mutate(Stat = "Total_offspring") %>% 
  rbind(mothers_bad_assignment %>% 
          group_by(Known_mother) %>% 
          summarize(count = n()) %>% 
          mutate(Stat = "Bad_assignment")) %>% 
  pivot_wider(names_from = "Stat", values_from = count) %>% 
  mutate(Bad_assignment = ifelse(is.na(Bad_assignment), 0, Bad_assignment),
         Prop_offspring_bad_assigned_per_brood = Bad_assignment / Total_offspring)

# Plot the relationship between the number of bad assignments and the brood size
relationship_brood_size_bad_attribution %>% 
  ggplot(aes(x = Total_offspring, y = Bad_assignment)) +
  geom_point()

# See if there is a linear relationship
lm(data = relationship_brood_size_bad_attribution, Bad_assignment ~ Total_offspring) %>% 
  summary()

# Look at the distribution of the proportion of bad assignments
relationship_brood_size_bad_attribution %>% 
  ggplot(aes(x = Prop_offspring_bad_assigned_per_brood)) +
  geom_histogram()

# Filter the mothers that have a high proportion (compared to their brood size)
# of bad assignments
mother_high_proportion_bad_assignments <- relationship_brood_size_bad_attribution %>% 
  filter(Prop_offspring_bad_assigned_per_brood > 0.5) %>% 
  pull(Known_mother)

maternity_comparison %>%
  filter(Known_mother %in% mother_high_proportion_bad_assignments)



########## Check replicates ##########
maternity_comparison %>% 
  filter(grepl("r", ID_DNA_RAD)) %>% 
  mutate(Sample_ID = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1],
         Replicate = ifelse(str_split_fixed(ID_DNA_RAD, "r", 2)[, 2] == 1, "replicate_1", "replicate_2")) %>% 
  select(-c(contains("2"), Is_mother_well_infered, Mother_ID, ProbMum1, ID_DNA_RAD, Known_mother)) %>% 
  pivot_wider(names_from = "Replicate", values_from = "InferredMum1") %>% 
  mutate(Is_mother_replicate_well_inferred = replicate_1 == replicate_2) %>% 
  filter(!Is_mother_replicate_well_inferred)
