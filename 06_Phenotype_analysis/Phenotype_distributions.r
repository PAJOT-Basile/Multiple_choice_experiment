##############################
########## Libraries #########
##############################
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "ggpubr", "ggh4x", "patchwork")
pacman::p_load(char = libraries, character.only = TRUE)
##############################
###### Useful variables ######
##############################
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")

my_theme <- theme_bw() +
  theme(text = element_text(size = 10))
##############################
###### Useful functions ######
##############################
get_number_indivs <- function(family_level, species, meta = metadata){
  meta %>% 
    filter(Sex == "M",
           Species == species,
           Family_level == family_level) %>% 
    nrow() %>% 
    return()
}

str_split_last <- function(string, pattern){
  nb_substrings <- str_count(string, pattern) + 1
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

"%not_in%" <- function(x, y){return(!(x %in% y))}

##############################
######### Import data ########
##############################
# Metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE)%>% 
  # Remove the dulicates
  filter(!grepl("_0", ID_DNA_RAD)) %>% 
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
  unique()


####################################
##### Total trait distribution #####
####################################
# Plot the distribution of the traits
nb_spines <- metadata %>% 
  select(Species, Sex, Family_level, P7.ep, P6.ep, P5.ep) %>% 
  filter(Sex == "M") %>%
  rowwise() %>% 
  pivot_longer(cols = starts_with("P"), values_to = "count", names_to = "Trait") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 10, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Species ~ Trait,
              labeller = as_labeller(c(
                "forsmani" = paste0("J. forsmani (n = ", get_number_indivs("parents", "forsmani"), " parents and n = ", get_number_indivs("offspring", "forsmani"), " offspring)"),
                "praehirsuta" = paste0("J. praehirsuta (n = ", get_number_indivs("parents", "praehirsuta"), " parents and n = ", get_number_indivs("offspring", "praehirsuta"), " offspring)"),
                "P5.ep" = "Pereiopod 5",
                "P6.ep" = "Pereiopod 6",
                "P7.ep" = "Pereiopod 7"))) +
  labs(x = "Number of spines",
       y = "Density",
       tag = "(B)") +
  theme_bw() +
  theme(text = element_text(size = 30),
        legend.position = "none",
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), units = "lines"))
# Plot the distribution of the traits
nb_setae <- metadata %>% 
  select(Species, Sex, Family_level, P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae) %>% 
  filter(Sex == "M") %>%
  rowwise() %>% 
  pivot_longer(cols = c(starts_with("P")), values_to = "count", names_to = "Trait") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Species ~ Trait,
              labeller = as_labeller(c(
                "forsmani" = paste0("J. forsmani (n = ", get_number_indivs("parents", "forsmani"), " parents and n = ", get_number_indivs("offspring", "forsmani"), " offspring)"),
                "praehirsuta" = paste0("J. praehirsuta (n = ", get_number_indivs("parents", "praehirsuta"), " parents and n = ", get_number_indivs("offspring", "praehirsuta"), " offspring)"),
                "P1.curv.setae" = "Pereiopod 1",
                "P2.curv.setae" = "Pereiopod 2",
                "P3.curv.setae" = "Pereiopod 3",
                "P4.curv.setae" = "Pereiopod 4")),
              scales = "free") +
  labs(x = "Number of setae",
       y = "Density",
       tag = "(A)") +
  theme_bw() +
  theme(text = element_text(size = 30),
        legend.position = "none",
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), units = "lines"))

size <- metadata %>%
  filter(Sex == "M") %>%
  select(Family_level, Size, Species) %>% 
  pivot_longer(Size, names_to = "Trait", values_to = "count") %>% 
  drop_na() %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Species ~ Trait,
              labeller = as_labeller(c(
                "Size" = "Size",
                "forsmani" = paste0("J.~forsmani~(n[parents]== ", get_number_indivs("parents", "forsmani"), "~and~n[offspring]==", get_number_indivs("offspring", "forsmani"), ")"),
                "praehirsuta" = paste0("J.~praehirsuta~(n[parents]== ", get_number_indivs("parents", "praehirsuta"), "~and~n[offspring]==", get_number_indivs("offspring", "praehirsuta"), ")")),
                default = label_parsed)) +
  labs(x = "Size (in mm)",
       y = "Density",
       tag = "(C)") +
  theme_bw() +
  theme(text = element_text(size = 30),
        legend.position = "none",
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 2, 0, 0), units = "lines"),
        axis.title.y = element_blank())



ggarrange(nb_setae,
  ggarrange(nb_spines, size,
            widths = c(3, 1)),
  ncol = 1) #%>% 
  # ggsave(plot = ., filename = "/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Figures/Phenotypic_traits/Distribution_traits.png",
  #        scale = 4, width = 2400, height = 1200, units = "px")

####################################
##### Distribution per segment #####
####################################
pivoted_data_males <- metadata %>% 
  filter(Sex == "M") %>%
  select(-c(Size, nb_tot_juv_brood1, nb_tot_juv_brood2, nb_juv_brood1,
            nb_juv_brood2, Name_offsprings_brood1, Name_offsprings_brood2,
            P7.ep, P7.ep.mid, P6.ep, P6.ep.dist, P5.ep, P5.ep.dist, P5.curv.setae,
            P4.curv.setae, P3.curv.setae, P2.curv.setae, P1.curv.setae,
            Mother_ID, Isolation_date, Death_date, Photo_ID, Note, Extracted)) %>% 
  pivot_longer(cols = -c(ID_DNA_RAD, Label, Alcool_ID, Sex, Phenotype, 
                         Family_level, Species),
               values_to = "count", names_to = "Trait") %>% 
  drop_na() %>%
  mutate(Pereiopod = str_split_fixed(Trait, "\\.", 2)[, 1],
         Segment = case_when(grepl("curv", Trait) ~ str_split_fixed(Trait, "\\.", 3)[, 2],
                             grepl("dist", Trait) ~ "Middle_carpus",
                             grepl("ep", Trait) ~ "carp",
                             TRUE ~ NA) %>% 
           factor(levels = c("isch", "merus", "Middle_carpus", "carp", "prop")),
         Setae_or_spine = ifelse(grepl("curv", Trait), "Setae", "Spines"),
         Size = case_when(grepl("larg", Trait) & !grepl("very", Trait) ~ "L",
                          grepl("very.larg", Trait) ~ "XL",
                          grepl("med", Trait) ~ "M",
                          grepl("small", Trait) ~ "S",
                          TRUE ~ NA) %>% 
           factor(levels = c("XL", "L", "M", "S")))

#### Look at the distribution of setae
### praehirsuta
# prae setae on P1-P4
pivoted_data_males %>% 
  filter(Setae_or_spine == "Setae",
         Species == "praehirsuta",
         Pereiopod %not_in% c("P6", "P5")) %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Segment ~ Pereiopod, scales = "free") +
  labs(x = "Number of setae") +
  my_theme

# prae setae on P5-P7
pivoted_data_males %>% 
  filter(Setae_or_spine == "Setae",
         Species == "praehirsuta",
         Pereiopod %in% c("P6", "P5")) %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Segment ~ Pereiopod, scales = "free") +
  labs(x = "Number of setae") +
  my_theme

### forsmani
# forsmani setae on P1-P4
pivoted_data_males %>% 
  filter(Setae_or_spine == "Setae",
         Species == "forsmani",
         Pereiopod %not_in% c("P6", "P5")) %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Segment ~ Pereiopod, scales = "free") +
  labs(x = "Number of setae") +
  my_theme

# forsmani setae on P5-P7
pivoted_data_males %>% 
  filter(Setae_or_spine == "Setae",
         Species == "forsmani",
         Pereiopod %in% c("P6", "P5")) %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Segment ~ Pereiopod, scales = "free") +
  labs(x = "Number of setae") +
  my_theme


#### Look at the distribution of spines
pivoted_data_males %>% 
  filter(Setae_or_spine == "Spines",
         Segment != "Middle_carpus") %>% 
  ggplot(aes(x = count, alpha = Family_level, fill = Species)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, position = "dodge") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
  facet_grid2(Segment ~ Pereiopod, scales = "free") +
  labs(x = "Number of setae") +
  my_theme
############################################
##### Look at variance in setae/spines #####
############################################
pivoted_data_males %>%
  group_by(Species, Pereiopod, Setae_or_spine, Segment, Family_level) %>%
  summarize(var = var(count)) %>%
  ggplot(aes(x = Segment, y = log10(var), fill = Species, alpha = Setae_or_spine)) +
  geom_col() +
  geom_hline(aes(yintercept = 0), colour = "black") +
  scale_fill_manual(values = colours_species) +
  scale_alpha_manual(values = c("Setae" = 0.5, "Spines" = 0.9)) +
  facet_grid2(Family_level ~ Pereiopod, scales = "free") +
  my_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

########################################
##### Correlation between segments #####
########################################
# As there is nearly no variation in the number of spines that are located on the middle of the carpus, we remove them from the analysis
pivoted_data_males <- pivoted_data_males %>% 
  filter(Segment != "Middle_carpus") %>% 
  mutate(Segment = case_when(Setae_or_spine == "Setae" ~ Segment,
                             TRUE ~ Size)) %>% 
  select(-Size) %>% 
  left_join(metadata %>% 
              select(ID_DNA_RAD, Size),
            by = "ID_DNA_RAD")

# Useful functions to plot
get_triangle_layout <- function(nb_rows){
  all_combinations <- get_all_combinations(1:nb_rows)
  text_layout <- ""
  row_nb <- 1
  nb_cols <- 1
  nb_empty_blocks <- nb_rows
  for (i in 1:nrow(all_combinations)){
    if (i == nb_cols){
      if (i == nrow(all_combinations)){
        text_layout <- paste0(text_layout, LETTERS[i])
        break
      }
      nb_empty_blocks <- nb_rows - row_nb
      empty_blocks <- paste0(rep("#", nb_empty_blocks), collapse = "")
      text_layout <- paste0(text_layout, LETTERS[i], empty_blocks, "\n")
      row_nb <- row_nb + 1
      nb_cols <- nb_cols + row_nb
    }else{
      text_layout <- paste0(text_layout, LETTERS[i])
    }
  }
  return(text_layout)
}
get_all_combinations <- function(x){
  all_combinations <- data.frame()
  for (i in seq_along(x)){
    for (j in i:length(x)){
      all_combinations <- add_table_to_df_in_iteration(all_combinations,
                                                       data.frame("x" = x[i], "y" = x[j]))
    }
  }
  return(all_combinations)
}
get_tag_letter <- function(x, y, vect){
  return(LETTERS[which(paste(x, y, sep = "_") == rev(vect))])
}
get_axis_label <- function(x){
  x_label <- case_when(
    x == "prop" ~ "Propode",
    x == "carp" ~ "Carpus",
    x == "isch" ~ "Ischium",
    x == "merus" ~ "Merus",
    x == "XL" ~ "Very large",
    x == "L" ~ "Large",
    x == "M" ~ "Medium",
    x == "S" ~ "Small",
    x == "Size" ~ "Size",
    TRUE ~ NA
  )
  return(x_label)
}
plot_correlation <- function(df, x, y, x_label, y_label, tag_letter, correlation_value){
  p <- df %>% 
    ggplot(aes(x = !!sym(x), y = !!sym(y), colour = Species, alpha = Family_level)) +
    geom_point(size = 4) +
    geom_smooth(se = FALSE, method = "lm") +
    annotate("text",
             x = min(df[[x]], na.rm = T) + 0.05 * max(df[[x]], na.rm = T),
             y = max(df[[y]], na.rm = T) - 0.1 * max(df[[y]], na.rm = T),
             label = correlation_value, size = 5) +
    scale_colour_manual(values = colours_species) +
    scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
    labs(x = x_label,
         y = y_label,
         tag = paste0("(", tag_letter, ")")) +
    my_theme
  return(p)
}
plot_histogram <- function(df, x, x_label, tag_letter){
  p <- df %>% 
    ggplot(aes(x = !!sym(x), fill = Species, alpha = Family_level)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20) +
    scale_fill_manual(values = colours_species) +
    scale_alpha_manual(values = c("parents" = 0.9, "offspring" = 0.3)) +
    labs(x = x_label,
         tag = paste0("(", tag_letter, ")")) +
    my_theme
  return(p)
}
get_blank_plot <- function(x_label, y_label, tag_letter){
  p <- ggplot() +
    geom_blank() +
    labs(x = x_label,
         y = y_label,
         tag = paste0("(", tag_letter, ")")) +
    my_theme
  return(p)
}
plot_pairs_of_segments <- function(df, x, y, vect){
  # Get the tag and axis names
  tag_letter <- get_tag_letter(y, x, vect)
  x_label <- get_axis_label(x)
  y_label <- get_axis_label(y)

    if (is.null((df[[x]])) | is.null(df[[y]])){
      p <- get_blank_plot(x_label, y_label, tag_letter)
    }else if (x == y){
      p <- plot_histogram(df, x, x_label, tag_letter)
    }else{
      # Remove NAs for the correlation
      df_no_na <- df %>% 
        drop_na()
      correlation_value <- cor(df_no_na[[x]], df_no_na[[y]]) %>% round(digits = 2)
      # Plot
      p <- plot_correlation(df, x, y, x_label, y_label, tag_letter, correlation_value)
      }
  return(p)
}

plot_correlations_pereiopod <- function(df, pereiopod, species, setae_or_spines){
  if (setae_or_spines %not_in% c("Setae", "Spines")){
    stop("Please write either 'Spines' or 'Setae' in the setae_or_spines argument")
  }
  # Filter the data
  df_temp <- df %>% 
    filter(Setae_or_spine == setae_or_spines,
           Species %in% species,
           Pereiopod %in% pereiopod)
  
  # Get all the combinations
  unique_segments <- df_temp %>% 
    pull(Segment) %>% 
    unique() %>% 
    as.character() %>% 
    c(., "Size")
  
  pairs_to_do <- get_all_combinations(unique_segments) %>% 
    mutate(Pairs = paste(x, y, sep = "_"))

  df_temp <- df_temp %>%
    select(-Trait) %>% 
    pivot_wider(names_from = "Segment", values_from = "count")
  
  # Make all the plots
  all_plots <- vector(mode = "list", length = nrow(pairs_to_do))
  for (row_nb in 1:nrow(pairs_to_do)){
    row_df_to_take <- nrow(pairs_to_do) - row_nb + 1
    all_plots[[row_nb]] <- plot_pairs_of_segments(df_temp, pairs_to_do[row_df_to_take, ]$y, pairs_to_do[row_df_to_take, ]$x, pairs_to_do$Pairs)
  }
  
  # Represent all the plots
  p <- all_plots[[1]]
  for (i in 2:length(all_plots)){
    p <- p + all_plots[[i]]
  }
  
  # Make the layout
  nb_rows <- length(unique_segments)
  
  layout_design <- get_triangle_layout(nb_rows)
  
  p <- p +
    plot_layout(design = layout_design, guides = "collect")
  return(p)
}

###### Number of setae
#### Number of setae prae P1
pivoted_data_males %>% 
  filter(Family_level == "parents") %>% 
  plot_correlations_pereiopod("P7", "forsmani", "Spines")