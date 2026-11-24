#' This script is written to look at the correlation of size with the other traits
#' studied here (the total number of spines and the total number of sexual setae)
#' 
#' Basile Pajot
#' Last updated: 18/11/2026


########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse", "ggExtra", "patchwork")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = libraries, character.only = TRUE)

########################################################
################### Useful variables ###################
########################################################
colours_species <- c("praehirsuta" = "navy", "forsmani" = "#3A9AB2")
colours_sexes <- c("Female" = "#D55E00", "Male" = "purple4")

my_theme <- theme_bw(base_family = "sans") +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


########################################################
##################### Import data ######################
########################################################
# Import the metadata
metadata <- read.table("Associated_data/metadata.tsv",
                       sep = "\t", header = TRUE) %>% 
  # Remove controls that
  filter(!grepl("_0", ID_DNA_RAD),
         Family_level == "parents") %>% 
  # remove sequencing duplicates
  mutate(ID_DNA_RAD = str_split_fixed(ID_DNA_RAD, "r", 2)[, 1] %>% 
           str_remove_all(., "a|b"),
         animal = as.integer(str_remove_all(ID_DNA_RAD, "E"))) %>% 
  unique() %>% 
  filter(!is.na(ID_DNA_RAD)) %>% 
  group_by(Sex, Species) %>% 
  mutate(across(c(Size, starts_with("Sum_")), ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE), .names = "St_{col}")) %>% 
  ungroup()


########################################################
###################### Plot data #######################
########################################################
max_size <- metadata %>% 
  filter(Sex == "M") %>% 
  pull(Size) %>% 
  max(na.rm = TRUE)

min_size <- metadata %>% 
  filter(Sex == "M") %>% 
  pull(Size) %>% 
  min(na.rm = TRUE)

fit_regression <- function(df, trait1, trait2){
  summ_model <- lm(formula(paste0(trait1, " ~ ", trait2)), data = df) %>% 
    summary()
  
  r_squared <- summ_model$r.squared
  beta_reg <- summ_model$coefficients[2, 1]
  intercept <- summ_model$coefficients[1, 1]
  p_val <- summ_model$coefficients[2, 4]
  
  tr1 <- case_when(trait1 == "Sum_curved_setaeP1P5_default" ~ "Setae",
                   trait1 == "Sum_spinesP6P7" ~ "Spines",
                   TRUE ~ "Size")
  tr2 <- case_when(trait2 == "Sum_curved_setaeP1P5_default" ~ "Setae",
                   trait2 == "Sum_spinesP6P7" ~ "Spines",
                   TRUE ~ "Size")
  
  if(intercept < 0){
    fmt <- "atop(R^2 == %.2f *','~ p == %.2e , %s == %.2f*' * '*%s~%.2f)"
  }else{
    fmt <- "atop(R^2 == %.2f *','~ p == %.2e , %s == %.2f*' * '*%s~+~%.2f)"
  }
  lab <- sprintf(fmt, r_squared, p_val, tr1, beta_reg, tr2, intercept)
  
  return(lab)
  
}

model_mpse <- metadata %>% 
  filter(Sex == "M", Species == "praehirsuta",
         !is.na(Size), !is.na(Sum_curved_setaeP1P5_default)) %>% 
  fit_regression("Size", "Sum_curved_setaeP1P5_default")

model_mpsp <- metadata %>% 
  filter(Sex == "M", Species == "praehirsuta",
         !is.na(Size), !is.na(Sum_spinesP6P7)) %>% 
  fit_regression("Size", "Sum_spinesP6P7")

model_mfse <- metadata %>% 
  filter(Sex == "M", Species == "forsmani",
         !is.na(Size), !is.na(Sum_curved_setaeP1P5_default)) %>% 
  fit_regression("Size", "Sum_curved_setaeP1P5_default")

model_mfsp <- metadata %>% 
  filter(Sex == "M", Species == "forsmani",
         !is.na(Size), !is.na(Sum_spinesP6P7)) %>% 
  fit_regression("Size", "Sum_spinesP6P7")


setae_f <- (metadata %>%
         filter(Sex == "M", Species == "forsmani") %>% 
         ggplot(aes(x = Size, y = Sum_curved_setaeP1P5_default, colour = Species)) +
         geom_point(size = 4, alpha = 0.7) +
         geom_smooth(method = "lm", se = FALSE, lwd = 1.1) +
           annotate("text", x = min_size, y = Inf, label = model_mfse,
                    size = 5, parse = TRUE, vjust = 1, hjust = 0) +
           scale_colour_manual(values = colours_species,
                             labels = c("forsmani" = "J. forsmani",
                                        "praehirsuta" = "J. praehirsuta"),
                             guide = "none") +
         xlim(c(min_size, max_size)) +
         my_theme +
         theme(legend.text = element_text(face = "italic"),
               axis.title = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank())) %>% 
  ggMarginal(fill = "#3A9AB2", type = "histogram", bins = 20)

setae_p <- (metadata %>%
         filter(Sex == "M", Species == "praehirsuta") %>% 
         ggplot(aes(x = Size, y = Sum_curved_setaeP1P5_default, colour = Species)) +
         geom_point(size = 4, alpha = 0.7) +
         geom_smooth(method = "lm", se = FALSE, lwd = 1.1) +
           annotate("text", x = max_size, y = Inf, label = model_mpse,
                    size = 5, parse = TRUE, vjust = 1, hjust = 1) +
           scale_colour_manual(values = colours_species,
                             labels = c("forsmani" = "J. forsmani",
                                        "praehirsuta" = "J. praehirsuta"),
                             guide = "none") +
         xlim(c(min_size, max_size)) +
         labs(x = "Size\n(mm)",
              y = "Number of curved setae") +
         my_theme +
         theme(legend.text = element_text(face = "italic"),
               axis.title = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank())) %>% 
  ggMarginal(fill = "navy", type = "histogram", bins = 20)


max_spines <- metadata %>% 
  filter(Sex == "M") %>% 
  pull(Sum_spinesP6P7) %>% 
  max(na.rm = TRUE)

min_spines <- metadata %>% 
  filter(Sex == "M") %>% 
  pull(Sum_spinesP6P7) %>% 
  min(na.rm = TRUE)


spines_f <- (metadata %>%
         filter(Sex == "M", Species == "forsmani") %>% 
         ggplot(aes(x = Size, y = Sum_spinesP6P7, colour = Species)) +
         geom_point(size = 4, alpha = 0.7) +
         geom_smooth(method = "lm", se = FALSE, lwd = 1.1) +
           annotate("text", x = max_size, y = min_spines, label = model_mfsp,
                    size = 5, parse = TRUE, vjust = 0, hjust = 1) +
           scale_colour_manual(values = colours_species,
                             labels = c("forsmani" = "J. forsmani",
                                        "praehirsuta" = "J. praehirsuta"),
                             guide = "none") +
         xlim(c(min_size, max_size)) +
         ylim(c(min_spines, max_spines)) +
         labs(x = "Size\n(mm)",
              y = "Number of spines") +
         my_theme +
         theme(legend.text = element_text(face = "italic"),
               axis.title.y = element_blank())) %>% 
  ggMarginal(fill = "#3A9AB2", type = "histogram", bins = 20)

spines_p <- (metadata %>%
         filter(Sex == "M", Species == "praehirsuta") %>% 
         ggplot(aes(x = Size, y = Sum_spinesP6P7, colour = Species)) +
         geom_point(size = 4, alpha = 0.7) +
         geom_smooth(method = "lm", se = FALSE, lwd = 1.1) +
           annotate("text", x = min_size, y = Inf, label = model_mpsp,
                    size = 5, parse = TRUE, vjust = 1, hjust = 0) +
           scale_colour_manual(values = colours_species,
                             labels = c("forsmani" = "J. forsmani",
                                        "praehirsuta" = "J. praehirsuta"),
                             guide = "none") +
         xlim(c(min_size, max_size)) +
         ylim(c(min_spines, max_spines)) +
         labs(x = "Size\n(mm)",
              y = "Number of spines") +
         my_theme +
         theme(legend.text = element_text(face = "italic"),
               axis.title.y = element_blank())) %>% 
  ggMarginal(fill = "navy", type = "histogram", bins = 20)


ann1 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~forsmani)"),
            parse = TRUE, size = 6) +
  theme_void()
ann2 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "italic(J.~praehirsuta)"),
            parse = TRUE, size = 6) +
  theme_void()

ann3 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Number~of~setae)"),
            parse = TRUE, size = 6, angle = 90) +
  theme_void()
ann4 <- ggplot() +
  geom_text(aes(x = 0, y = 0, label = "bold(Number~of~spines)"),
            parse = TRUE, size = 6, angle = 90) +
  theme_void()




size_trait_relationship <- ann1 + ann2 +
  ann3 + setae_f + setae_p +
  ann4 + spines_f + spines_p +
  plot_layout(design = "#AB
              CDE
              FGH",
              widths = c(0.25, 1, 1),
              heights = c(0.15, 0.8, 1))

size_trait_relationship
ggsave(plot = size_trait_relationship,
       filename = "../Figures/Size_trait_covariance.png",
       width = 800, height = 800, units = "px", scale = 4)


########################################################
################# Regression on  data ##################
########################################################
# Size and number of setae:
metadata %>% 
  filter(Sex == "M", Species == "praehirsuta",
         !is.na(Size)) %>% 
  lm(data = ., Sum_curved_setaeP1P5_default ~ Size) %>% 
  summary()
