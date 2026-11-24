#' This script is written to look at the correlation between the number of
#' spines that was counted with and without the carps. As the carps are not always
#' counted with a large acuracy, we remove them from the total number of setae 
#' and see if there is a correlation.
#' 
#' Basile Pajot
#' Last updated: 17/11/2026


########################################################
####################### Libraries ######################
########################################################
libraries <- c("tidyverse")
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
######################## Model #########################
########################################################
summ_model <- metadata %>% 
  filter(Species == "praehirsuta", Sex == "M") %>% 
  mutate(Sum_curved_setae_without_carp = Sum_curved_setaeP1P5_default - default_P5.carp.curv.setae - default_P4.carp.curv.setae - default_P3.carp.curv.setae - default_P2.carp.curv.setae - default_P1.carp.curv.setae) %>% 
  filter(!is.na(Sum_curved_setaeP1P5_default)) %>% 
  lm(data = ., Sum_curved_setaeP1P5_default ~ Sum_curved_setae_without_carp) %>% 
  summary()

r_squared <- summ_model$r.squared
beta_reg <- summ_model$coefficients[2, 1]
intercept <- summ_model$coefficients[1, 1]
p_val <- summ_model$coefficients[2, 4]

if(intercept < 0){
  fmt <- "atop(R^2 == %.2f *','~ p == %.2e , Sum_with_carp == %.2f*' * '*Sum_without_carp~%.2f)"
}else{
  fmt <- "atop(R^2 == %.2f *','~ p == %.2e , Sum_with_carp == %.2f*' * '*Sum_without_carp~+~%.2f)"
}
lab <- sprintf(fmt, r_squared, p_val, beta_reg, intercept)
########################################################
###################### Plot data #######################
########################################################
relationship_with_without_carp <- metadata %>% 
  filter(Species == "praehirsuta", Sex == "M") %>% 
  mutate(Sum_curved_setae_without_carp = Sum_curved_setaeP1P5_default - default_P5.carp.curv.setae - default_P4.carp.curv.setae - default_P3.carp.curv.setae - default_P2.carp.curv.setae - default_P1.carp.curv.setae) %>% 
  ggplot(aes(x = Sum_curved_setae_without_carp, y = Sum_curved_setaeP1P5_default, colour = Species)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(se = FALSE, method = "lm") +
  annotate("text", x = Inf, y = Inf, label = lab,
           hjust = 1.2, vjust = 1.5, size = 7, parse = TRUE) +
  scale_colour_manual(values = colours_species,
                      labels = c("forsmani" = "J. forsmani",
                                 "praehirsuta" = "J. praehirsuta"),
                      guide = "none") +
  labs(x = "Number of curved setae without carps",
       y = "Number of curved setae with carps") +
  facet_wrap(vars(Species), nrow = 1, scales = "free",
             labeller = as_labeller(c("forsmani" = "J. forsmani",
                                      "praehirsuta" = "J. praehirsuta"))) +
  scale_x_continuous(limits = c(0, NA)) +
  my_theme +
  theme(legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic"))

relationship_with_without_carp
ggsave(filename = "../Figures/Suppfigs/Setae_with_without_carps.png",
       relationship_with_without_carp,
       scale = 4, width = 1000, height = 600, units = "px")
