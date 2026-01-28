##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet")
pacman::p_load(char = libraries, character.only = TRUE)

##############################
########### Useful ###########
##############################
my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

# Function not in
"%!in%" <- function(x, y){ return(!(x %in% y))}

##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  mutate(Size = str_replace_all(Size, ",", ".") %>% as.numeric)

parents <- metadata %>% 
  filter(Family_level == "parents") %>% 
  pull(ID_DNA_RAD)


# Import relatedness between all individuals
relatedness_formsani_females <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/08_Calling/ddRAD_multiple_choice_experiment/relatedness/Relatedness_female_parents_forsmani.tsv",
                                   sep = "\t", header = TRUE)


relatedness_formsani_females %>%
  ggplot(aes(x = RELATEDNESS_PHI)) + 
  geom_histogram(bins = 100)


relatedness_formsani_females %>% 
  filter(RELATEDNESS_PHI > 0.1 & RELATEDNESS_PHI < 0.4)



#################################"


relatedness_everyone %>% 
  filter(INDV1 %in% parents & INDV2 %in% parents) %>% 
  ggplot(aes(x = RELATEDNESS_AJK)) +
  geom_histogram(bins = 100)
# Find parents that are closely-related:
relatedness_everyone %>% 
filter(INDV1 %in% parents & INDV2 %in% parents,
       RELATEDNESS_AJK > 0.07) %>% 
mutate(identity = INDV1 == INDV2) %>% 
filter(!identity) %>%
left_join(metadata %>%
            select(ID_DNA_RAD, Label, PCA_species),
          by = join_by("INDV1" == "ID_DNA_RAD")) %>% 
  left_join(metadata%>%
              select(ID_DNA_RAD, Label, PCA_species),
            by = join_by("INDV2" == "ID_DNA_RAD"),
            suffix = c("_parent1", "_parent2")) %>% 
  arrange(PCA_species_parent1) %>% 
  group_by(PCA_species_parent1) %>% 
  summarize(count = n())
