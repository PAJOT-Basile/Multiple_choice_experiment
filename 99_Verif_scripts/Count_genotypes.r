########## Libraries ##########
if (!require("pacman")) install.packages("pacman")
libraries <- c("tidyverse", "vcfR", "adegenet")
pacman::p_load(char = libraries, character.only = TRUE)

########### Useful ###########
geom_manhattan <- function(df, mapping, thresholding = FALSE, absolute = TRUE, palette = c("grey41", "grey4"), filter_low_freqs = TRUE, ...){
  #' Draw manhattan plot
  #' 
  #' This function traces a manhattan plot relying on the ggplot2 aesthetics.
  #' Some additional arguments have been added to be able to personalise the plots
  #' to the best.
  #' 
  #'@param df (data.frame).
  #'    This data frame contains the information to plot. It requires at least
  #'    two columns: the position in the genome and the column to plot along the
  #'    genome.
  #'@param mapping (mapping object).
  #'    This is the mapping of the graph just like in the ggplot2 library. The
  #'    mapping in this function does not require an x column because the position
  #'    along the genome is taken by default.
  #'@param thresholding (boolean). (default = FALSE)
  #'    This is an argument that changes the output of the function. It outputs
  #'    the graph and a maximum cumulative position along the whole genome. It is
  #'    to be used with the "thresholds_manhattan" function.
  #'@param absolute (boolean). (default = TRUE)
  #'    This arguments simply says if you want to plot the absolute values in the
  #'    selected column or if you plot the real value
  #'@param palette (vector). (default = c("grey71", "orange2"))
  #'    This argument is the colour palette to use to distinguish between the 
  #'    successive chromosomes. 
  #'@param ...
  #'    In these arguments, you can add any arguments that you would give a 
  #'    ggplot2 graph outside of the aesthetics.
  #'    
  #'@returns (ggplot2 object)
  #'    This returns the manhattan plot of the required column values along the
  #'    genome
  #'@returns (numeric value)
  #'    If the "thresholding" argument is TRUE, then this function also returns
  #'    the maximum cumulative position in the whole genome.
  #'@section Warning:
  #' This function has been made to deal only with two-coloured palettes.
  #' @export
  
  # First, we added some error checking to keep the function from running for
  # nothing
  if ("Position" %!in% names(df)){
    stop("This function needs a column called 'Position' that contains the position of each SNP on the genome")
  }
  if (is.numeric(df$Position[1])){
    if ("Chromosome" %!in% names(df)){
      stop("If the positions are already usable, please name a column 'Chromosome' to be used as a chromosome reference.")
    }
  }
  if ("y" %!in% names(mapping)){
    stop("This function requires a y aesthetic.")
  }
  
  # We store the original mapping because the mapping will be modified in this 
  # function, so we need to be able to compare the modified to the original one
  mapping_ori <- mapping
  
  # Then, we extract the name of the column that we have to represent along the
  # y axis.
  colName <- mapping$y[2] %>% as.character
  # And we transform it into a name so that we can use it with the tidyverse.
  colName <- as.name(substitute(colName))
  
  # Now that all the argument importation has been done, we are going to select
  # the columns in the input dataframe that we want to keep in the analysis. To 
  # do this, we make a list of all the parameters that were called to use in the 
  # "matches" function of select
  # We initialise the string "to_select" to an empty character that will be
  # complementeds
  to_select <- ""
  # We iterate over the names of the arguments passed to the function (except
  # for the first one which is the name of the datafame)
  for (i in 1:length(mapping)){
    # We differentiate the first one because there will be no "|" character
    # before it
    if (to_select == ""){
      to_select <- quo_name(mapping[[i]])
    }else{
      # Add the name of the parameters from the mapping to keep
      to_select <- paste0(to_select, "|", quo_name(mapping[[i]]))
    }
  }
  
  # Then, we separate cases where the "Position" column is already numeric
  # (position on one chromosome for example) from the case where the "Position"
  # column is under the adegenet format.
  if (is.numeric(df$Position[1])){
    # Now, we keep only the columns of interest in the dataframe
    To_plot <- df %>% 
      # We keep the columns of interest in the data frame
      select(Position, Chromosome, matches(to_select)) %>%
      # We re-order the levels of the "Chromosome" variable to have them in 
      # increasing order (from 1 to n)
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>%  
                        mutate(second = case_when(
                          grepl("\\.", Chromosome) ~ str_split_fixed(Chromosome, "\\.", 2)[, 2],
                          TRUE ~ letters[1]
                          )) %>%
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome)), second) %>%
                        pull(Chromosome) %>%
                        unique())) %>% 
      # We remove missing values if there are any
      drop_na()
  }else{
    # If the "Position" column is not numeric, it is likely they are at the 
    # adegent format
    To_plot <- df %>%
      # So, we use the created function to separate the position on the
      # chromosome and the chromosome name
      transform_position_ade2tidy() %>% 
      # Then, we re-order the levels of the "Chromosome" variable to have them in 
      # increasing order (from 1 to n)
      mutate(Chromosome = Chromosome %>% 
               factor(levels = df %>%
                        transform_position_ade2tidy() %>%
                        mutate(second = case_when(
                          grepl("\\.", Chromosome) ~ str_split_fixed(Chromosome, "\\.", 2)[, 2],
                          TRUE ~ letters[1]
                        )) %>%
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome)), second) %>%
                        pull(Chromosome) %>% 
                        unique() )) %>% 
      # And we select the columns of interest
      select(Position, Chromosome, matches(to_select)) %>% 
      # Finally, we drop the missing values
      drop_na()
  }
  
  # Then, we make a cumulative data frame with the positions of each end and beginning of chromosome
  data_cum <- To_plot %>% 
    group_by(Chromosome) %>% 
    mutate(second = case_when(
      grepl("\\.", Chromosome) ~ str_split_fixed(Chromosome, "\\.", 2)[, 2],
      TRUE ~ letters[1]
    )) %>%
    arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome)), second) %>%
    select(-second) %>% 
    summarise(max_bp = Position %>% max) %>% 
    # Here, we make a new column that contains the maximum position of the
    # previous chromosome to add this value to the position of the positions of
    # said chromosome
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(Chromosome, bp_add)
  
  # We join it to the filtered dataframe
  To_plot_manhat <- To_plot %>% 
    inner_join(data_cum, by="Chromosome") %>%
    # And we add the cumulative position for each chromosome along the genome
    # to place each chromosome one after the other
    mutate(bp_cum = Position + bp_add)
  
  # Then, we find the center of the chromosomes
  chromosome_centers <- To_plot_manhat %>% 
    group_by(Chromosome) %>%
    # The values of the centers of the chromosomes are approximated using the
    # mean function
    summarise(center = bp_cum %>% mean)
  
  
  # We select two colours that will be used to distinguish chromosomes
  colour_chromosome <- rep(palette, To_plot$Chromosome %>% unique %>% length)
  
  # Once all this is done, we create a vector of boolean to use as indicator to
  # take or not the absolute value of the column to plot
  absolute_list <- rep(absolute, nrow(To_plot_manhat))
  
  # Finally, we mutate, if needed, the column to plot with the absolute value
  To_plot_manhat <- To_plot_manhat  %>% 
    right_join(chromosome_centers, by="Chromosome") %>% 
    mutate(!!as.symbol(colName) := ifelse(absolute_list == TRUE, !!as.symbol(colName) %>% abs, !!as.symbol(colName)))
  if (filter_low_freqs){
    # We also filter some values to lighten the plot a little bit
    To_plot_manhat <- To_plot_manhat %>% 
      filter(!!as.symbol(colName) %>% abs > 0.05)
  }
  
  # Once this is done, we create the architecture of the plot (x axis, name of
  # the axis and the theme to use)
  p <- ggplot(data = To_plot_manhat, aes(x = bp_cum)) +
    # The x scale we use is just to differenciate the chromosomes (indicated by
    # their number)
    scale_x_continuous(labels = ((chromosome_centers$Chromosome %>%
                                    str_split_fixed(., "_", 2))[, 1] %>%
                                   str_split_fixed(., "G", 2))[, 2],
                       breaks = chromosome_centers$center) +
    labs(x = "Chromosomes",
         y = colName) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 20))
  
  # We are going to build a string containing the function to call in the ggplot.
  # To do this, we are going to make a list of the default values for aesthetic
  # parameters and we iterate over them.
  list_default_parameters <- list("y" = ".", "colour" = "Chromosome", "size" = 2)
  # We initialise the function to the geom_point function with the mapping we are
  # going to use
  function_to_call <- "geom_point(aes(mapping)"
  for (param in names(list_default_parameters)){
    default_param <- list_default_parameters[[param]]
    # For each parameter, if it is not in the mapping, we add it to the function
    # to call 
    if (param %!in% names(mapping)){
      if (param != "colour"){
        function_to_call <- paste0(function_to_call, ", ", param, " = ", default_param)
      }else{
        # The only exception to this is that if there is no colour mapping, we
        # add one to distinguish the chromosomes
        mapping <- c(mapping, aes(colour = Chromosome %>% as.factor))
        class(mapping) <- "uneval"
      }
    }else{
      # If the parameters are in the mapping, we add their names to the labels
      # of the plot
      p$labels[[param]] <- mapping[[param]] %>% quo_name
    }
  }
  
  # If there are some supplementary ggplot arguments they are added to the plot 
  if (length(list(...))){
    function_to_call <- paste0(function_to_call, ", ...")
  }
  # The function to call variable is completed and closed here.
  function_to_call <- paste0(function_to_call, ", inherit.aes = TRUE)")
  if ("colour" %in% names(mapping_ori)){
    # If there is a colour mapping in the plot, we can not distinguish
    # chromosomes by colouring the points of the plot by chromosome, so we create
    # polygons in the background that will be coloured differently.
    # Fisrt, we isolate the name of the colour column
    colour <- mapping$colour[2] %>% as.character
    # And we trasform it into a name to use it in the tidyverse
    colour <- as.name(substitute(colour))
    
    # We use the geom_box_background function to get the data frame of polygon
    # delimitations
    boxes <- geom_box_background(To_plot_manhat, colName, chromosome_centers)
    
    # We add the polygons to the graph
    p <- p +
      # First, the polygons so they are in the background
      geom_polygon(data = boxes, aes(x = x_value, y = y_value, fill = Chromosome), colour = NA) +
      # Then, we chose how to colour the chromosomes
      scale_fill_manual(values = colour_chromosome, guide="none")
    
    # We transform the beginning of the funnction to call variable so that it does
    # not inherit the table boxes as argument
    function_to_call <- function_to_call %>% str_replace("geom_point\\(", "geom_point\\(data = To_plot_manhat, ")
    
    # Then, we use a continuous colour scale if we have a continuous variable and a 
    # discrete colour scale if we have discrete values.
    if (is.continuous(To_plot_manhat[[as.character(colour)]])){
      p <- p +
        # We add the function to call to the graph
        eval(parse(text = function_to_call)) +
        # We use a continuous colour scale
        scale_colour_gradientn(name = as.character(colour),
                               colours = c("darkorchid4", "darkorchid", "mediumorchid1", "magenta"))
    }else{
      # In the case where the colour column is discrete, we have to separate the 
      # case where the colour column is the same as the one we want to plot the 
      # values of along the genome form the case where they are different.
      if (colour == colName){
        # If the two columns are the same, we have to separate them to transform
        # the colour column into factors. So, we create a duplicate column
        To_plot_manhat <- To_plot_manhat %>% 
          mutate(colour_column = !!as.symbol(colour) %>% factor(levels = df %>% 
                                                                  select(as.symbol(colour)) %>% 
                                                                  unique() %>% arrange(!!as.symbol(colour)) %>% 
                                                                  as.vector %>% unname %>% unlist))
        # In this case, we have to modify the mapping of the graph
        mapping$colour <- quo(colour_column)
        # And select some colours to use
        colour_polygon_chrom <- To_plot_manhat %>% 
          select(colour_column) %>%
          n_distinct() %>% 
          plasma()
      }else{
        # If the columns are different, the process is the same, except we simply
        # use the colour column rather than creating a duplicate column.
        To_plot_manhat <- To_plot_manhat %>% 
          mutate(!!as.symbol(colour) := !!as.symbol(colour) %>% factor(levels = df %>% 
                                                                         select(as.symbol(colour)) %>% 
                                                                         unique() %>% arrange(!!as.symbol(colour)) %>% 
                                                                         as.vector %>% unname %>% unlist))
        
        colour_polygon_chrom <- To_plot_manhat %>% 
          select(!!as.symbol(colour)) %>%
          n_distinct() %>% 
          plasma()
      }
      # We add this to the plot and change the colour scale
      p <- p +
        eval(parse(text = function_to_call)) +
        scale_colour_manual(name = as.character(colour),
                            values = colour_polygon_chrom,
                            drop = FALSE)
    }
    # Finally, we transform the layers of the plot so it can add the required 
    # mapping
    p$layers[[2]]$computed_mapping <- NULL
    p$layers[[2]]$mapping <- mapping
    
    
  }else{
    # If the colour argument is not in the mapping, we simply colour using the 
    # chromosomes (one colour per chromosome)
    p <- p +
      # We add the points to the plot
      eval(parse(text = function_to_call)) +
      # We change the colour scale
      scale_colour_manual(values = colour_chromosome, guide="none")
    
    # And in the same way, we modify the layers of the plot so as to use the 
    # required aesthetics
    p$layers[[1]]$computed_mapping <- NULL
    p$layers[[1]]$mapping <- mapping
    
  }
  # Finally, if there is a thresholding (i.e. if the function is embedded in the
  # thresholding function), then, we simply return a list of parameters
  if (thresholding) {
    return(list("plot" = p, "max_value" = To_plot_manhat %>% select(bp_cum) %>% max))
  }else{
    # Otherwise, we return the plot
    return(p)
  }
}

my_theme <- theme_bw() +
  theme(text = element_text(size = 20))

colours_species <- c("praehirsuta" = "navy",
                     "forsmani" = "#3A9AB2",
                     "hybrid" = "red",
                     "hybrid_genotype" = "green")

# Function not in
"%!in%" <- function(x, y){ return(!(x %in% y))}

Get_chromosome_summary_information <- function(chromosome_recap){
  #' This function makes a recap of the chromosome information
  #' (max position, center, ...).
  #' /!\ THE CHROMOSOMES HAVE TO BE IN THE CORRECT ORDER AND THEIR COLUMN HAS TO BE NAMED "Chromosome"
  if ("Chromosome" %!in% names(chromosome_recap)) stop("Is your column containing the names of the chromosomes named 'Chromosome'?")
  if ("Position" %!in% names(chromosome_recap)) stop("Is your column containing the positions along the chromosomes named 'Position'?")
  
  
  chromosome_recap %>% 
    group_by(Chromosome) %>% 
    summarize(max_position = max(Position, na.rm = TRUE)) %>% 
    mutate(Start_chromosome = cumsum(lag(max_position + 1, default = 0)),
           Center_chromosome = round(max_position) / 2 + Start_chromosome) %>% 
    return()
}

prepare_data_for_manhattan <- function(df, sum_stats = summary_stats, chromosome_info = chromosome_information, chrom_order = chromosome_order){
  df %>% 
    separate(loci, c("Locus", "Col", "direction_RAD"), sep = ":") %>% 
    select(-direction_RAD) %>% 
    mutate(Locus = Locus %>% as.numeric,
           Col = as.numeric(Col)) %>% 
    left_join(sum_stats %>% 
                select(Locus, Col, BP, Chrom) %>% 
                dplyr::rename(Chromosome = Chrom,
                              Position = BP),
              by = c("Locus", "Col")) %>% 
    left_join(chromosome_info, by = "Chromosome") %>% 
    mutate(Chromosome = Chromosome %>% 
             factor(levels = chrom_order)) %>% 
    return()
}

get_marker_types <- function(genetic_data, mom, dad, offspring){
  indivs_of_interest <- paste(c(mom, dad, offspring), collapse = "|")
  Type_markers <- genetic_data[grepl(indivs_of_interest, indNames(genetic_data))]@tab %>% 
    t() %>% 
    as_tibble() %>% 
    mutate(Position = colnames(genetic_data@tab)) %>% 
    drop_na() %>% 
    filter(grepl("\\.1", Position)) %>% 
    mutate(Mom_ID = mom,
           Dad_ID = dad,
           Offspring_ID = offspring)
  
  if (names(Type_markers)[2] != mom){
    Type_markers <- Type_markers %>% 
      rename(!!sym(mom) := names(Type_markers)[2])
  }
  
  Type_markers <- Type_markers %>% 
    mutate(Type = case_when(
      !!sym(mom) == !!sym(offspring) & !!sym(dad) == !!sym(offspring) ~ "Uninformative",
      !!sym(mom) == 1 & !!sym(dad) == 1 & !!sym(offspring) != 1 ~ "Uninformative",
      !!sym(mom) == 0 & !!sym(dad) == 0 & !!sym(offspring) != 0 ~ "Creationist",
      !!sym(mom) == 0 & !!sym(dad) == 1 & !!sym(offspring) != 2 ~ "Uninformative",
      !!sym(mom) == 0 & !!sym(dad) == 1 & !!sym(offspring) == 2 ~ "Uniparental_dad",
      !!sym(mom) == 0 & !!sym(dad) == 2 & !!sym(offspring) == 1 ~ "Uninformative",
      !!sym(mom) == 0 & !!sym(dad) == 2 & !!sym(offspring) == 0 ~ "Uniparental_mom",
      !!sym(mom) == 0 & !!sym(dad) == 2 & !!sym(offspring) == 2 ~ "Uniparental_dad",
      !!sym(mom) == 1 & !!sym(dad) == 0 & !!sym(offspring) != 2 ~ "Uninformative",
      !!sym(mom) == 1 & !!sym(dad) == 0 & !!sym(offspring) == 2 ~ "Uniparental_mom",
      !!sym(mom) == 1 & !!sym(dad) == 2 & !!sym(offspring) != 0 ~ "Uninformative",
      !!sym(mom) == 1 & !!sym(dad) == 2 & !!sym(offspring) == 0 ~ "Uniparental_mom",
      !!sym(mom) == 2 & !!sym(dad) == 0 & !!sym(offspring) == 1 ~ "Uninformative",
      !!sym(mom) == 2 & !!sym(dad) == 0 & !!sym(offspring) == 0 ~ "Uniparental_dad",
      !!sym(mom) == 2 & !!sym(dad) == 0 & !!sym(offspring) == 2 ~ "Uniparental_mom",
      !!sym(mom) == 2 & !!sym(dad) == 1 & !!sym(offspring) != 0 ~ "Uninformative",
      !!sym(mom) == 2 & !!sym(dad) == 1 & !!sym(offspring) == 0 ~ "Uniparental_dad",
      !!sym(mom) == 2 & !!sym(dad) == 2 & !!sym(offspring) != 2 ~ "Creationist",
      TRUE ~ NA
    )) %>% 
    rename(Mom = mom,
           Dad = dad,
           Offspring = offspring)
  
  Summary_table <- Type_markers %>% 
    group_by(Mom, Dad, Offspring) %>% 
    summarize(nb_SNPs = n(), .groups = "drop_last") %>%
    left_join(Type_markers %>% 
                group_by(Mom, Dad) %>% 
                summarize(total_snps = n(), .groups = "drop_last"),
              by = c("Mom", "Dad")) %>% 
    mutate(Proportion = nb_SNPs / total_snps,
           Expected_prop_Mendel = case_when(
             Mom == 0 & Dad == 0 & Offspring == 0 ~ 1,
             Mom == 0 & Dad == 1 & Offspring != 2 ~ 0.5,
             Mom == 0 & Dad == 2 & Offspring == 1 ~ 1,
             Mom == 1 & Dad == 0 & Offspring != 2 ~ 0.5,
             Mom == 1 & Dad == 1 & Offspring != 1 ~ 0.25,
             Mom == 1 & Dad == 1 & Offspring == 1 ~ 0.5,
             Mom == 1 & Dad == 2 & Offspring != 0 ~ 0.5,
             Mom == 2 & Dad == 0 & Offspring == 1 ~ 1,
             Mom == 2 & Dad == 1 & Offspring != 0 ~ 0.5,
             Mom == 2 & Dad == 2 & Offspring == 2 ~ 1,
             TRUE ~ 0
           ),
           Expected_nb_SNPs = round(Expected_prop_Mendel * total_snps)) %>% 
    relocate(Mom, Dad, Offspring, nb_SNPs, Expected_nb_SNPs, total_snps, Proportion, Expected_prop_Mendel)
  return(list(
    "Type_markers" = Type_markers,
    "Summary_table" = Summary_table
  ))
}

add_markers_to_manhattan <- function(type_markers, Fst, sumstats = summary_stats, chrom_info = chromosome_information, chrom_order = chromosome_order){
  markers_of_interest <- type_markers %>% 
    separate_wider_delim(Position, names = c("Locus", "Col", "Allele"), ":") %>% 
    mutate(across(c(Locus, Col), ~ as.numeric(.))) %>% 
    left_join(sumstats, by = c("Locus", "Col")) %>% 
    left_join(chrom_info, by = join_by("Chrom" == "Chromosome")) %>% 
    select(-Chromosome) %>%
    mutate(BP_cumul2 = BP + Start_chromosome) %>%
    left_join(Fst %>% 
                prepare_data_for_manhattan(sumstats, chrom_info, chrom_order),
              by = c("Col", "Locus", "max_position", "Start_chromosome", "Center_chromosome"))
  
  p <- (Fst %>% 
          prepare_data_for_manhattan(sumstats, chrom_info, chrom_order) %>% 
          geom_manhattan(aes(y = value), filter_low_freqs = FALSE)) +
    geom_point(data = markers_of_interest,
               aes(x = BP_cumul2, y = value), color = "red", alpha = 0.5) +
    labs(y = "Fst")
  return(p)
}
########### Import data ###########
# Import the metadata
metadata <- read.table("/shared/projects/sexisol/input/Basile/Multiple_choice_experiment/Data/metadata.tsv",
                       sep = "\t", header = TRUE)

# Thinned vcf
data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_whithout_unwanted_snps_no_error_and_good_indivs.vcf") %>% 
  vcfR2genind()

# Get the sizes of the chromosomes and compute where the next chromosome will start on the manhattan plot
chromosome_sizes <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/chromosome_sizes.tsv",
                               sep = "\t", header = TRUE)

# Import the sumstats of the thinned vcf (locus, chromosome, position, ...)
summary_stats <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/sumstats_thinned_50k.tsv",
                            sep = "\t", header = TRUE)

# The chromosome information needed to do a manhattan plot
chromosome_information <- summary_stats %>%
  select(Chrom, BP) %>%
  arrange(Chrom) %>%
  rename(Chromosome = Chrom, Position = BP) %>%
  Get_chromosome_summary_information()

# Import the Fst
load("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/Pop_gen/Fst_prae_forsmani.Rdata")

# Import the names of the hybrid offspring
load("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/tmp_colony/Hybrid_offsrping.Rdata")

########### Test for one hybrid offspring ###########
comb <- hybrid_offspring[3, ]
"E4731"
Type_markers <- get_marker_types(data, "E4914", "E4754", "E5258")

# Random individuals
get_marker_types(data, "E4914", "E4754", "E5258")$Type_markers %>% 
  group_by(Type) %>%
  summarize(count = n()) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~ "Total")))

# get_marker_types(data, "E4914", "E4754", "E5258")$Summary_table %>% 
#   ungroup() %>% 
#   mutate(across(c(Mom, Dad), ~ paste(., Offspring, sep = "_")),
#          across(c(Mom, Dad), ~ factor(., levels = c("0_0", "0_1", "0_2",
#                                                     "1_0", "1_1", "1_2",
#                                                     "2_0", "2_1", "2_2")))) %>%
#   pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
#   ggplot(aes(x = Mom, y = Dad, fill = Proportion)) +
#   geom_tile() +
#   facet_wrap(vars(Stat)) +
#   my_theme

# Intra-specific trio
get_marker_types(data, "E4919", "E4764", "E5245")$Type_markers %>% 
  group_by(Type) %>%
  summarize(count = n()) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~ "Total")))

# get_marker_types(data, "E4919", "E4764", "E5245")$Summary_table %>% 
#   ungroup() %>% 
#   mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
#          across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
#                                                     "1_0", "1_1", "1_2",
#                                                     "2_0", "2_1", "2_2")))) %>%
#   pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
#   ggplot(aes(x = Mom, y = Dad, fill = Proportion)) +
#   geom_tile() +
#   facet_wrap(vars(Stat)) +
#   my_theme

# Hybrid trio
get_marker_types(data, comb$Mom, comb$Dad, comb$Offspring)$Type_markers %>% 
  group_by(Type) %>%
  summarize(count = n()) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~ "Total")))

# get_marker_types(data, comb$Mom, comb$Dad, comb$Offspring)$Summary_table %>% 
#   ungroup() %>% 
#   mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
#          across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
#                                                     "1_0", "1_1", "1_2",
#                                                     "2_0", "2_1", "2_2")))) %>%
#   pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
#   ggplot(aes(x = Mom, y = Dad, fill = Proportion)) +
#   geom_tile() +
#   facet_wrap(vars(Stat)) +
#   my_theme


data <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_no_error_10_DP_100.vcf") %>% 
  vcfR2genind()
# get_marker_types(data, "E4919", "E4764", "E5245")$Summary_table %>% 
get_marker_types(data, comb$Mom, "E4749", "E5454")$Summary_table %>% 
  ungroup() %>% 
  # mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
  #        across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
  #                                              "1_0", "1_1", "1_2",
  #                                              "2_0", "2_1", "2_2")))) %>%
  pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
  filter(Stat == "Expected_prop_Mendel") %>% 
  mutate(Stat = "Proportion") %>% 
  mutate(Cross = "Expected") %>% 
  rbind(get_marker_types(data, "E4935", comb$Dad, comb$Offspring)$Summary_table %>%
          ungroup() %>%
          # mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
          #        across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
          #                                              "1_0", "1_1", "1_2",
          #                                              "2_0", "2_1", "2_2")))) %>%
          pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>%
          mutate(Cross = "Heterospecific")) %>% 
  rbind(get_marker_types(data, "E4894", "E4733", "E5017")$Summary_table %>% 
          ungroup() %>% 
            # mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
            #        across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
            #                                                   "1_0", "1_1", "1_2",
            #                                                   "2_0", "2_1", "2_2")))) %>%
          pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
          mutate(Cross = "Forsmani")) %>% 
  rbind(get_marker_types(data, "E4907b", "E4748", "E5142")$Summary_table %>% 
          ungroup() %>% 
          # mutate(across(c(Mom), ~ paste(., Offspring, sep = "_")),
          #        across(c(Mom), ~ factor(., levels = c("0_0", "0_1", "0_2",
          #                                                   "1_0", "1_1", "1_2",
          #                                                   "2_0", "2_1", "2_2")))) %>%
          pivot_longer(cols = contains("rop"), names_to = "Stat", values_to = "Proportion") %>% 
          mutate(Cross = "Praehirsuta")) %>% 
  filter(Stat == "Proportion") %>%
  ggplot() +
  geom_tile(aes(x = Mom, y = Dad, fill = Proportion)) +
  scale_fill_gradient(low = "gold", high ="orangered1") + 
  ggh4x::facet_grid2(Cross ~ Offspring) +
  my_theme

get_marker_types(data, "E4894", "E4733", "E5017")$Summary_table %>% view()

get_marker_types(data, "E4894", "E4733", "E5017")$Type_markers %>% 
  mutate(Position = str_remove_all(Position, "\\.1")) %>% 
  left_join(error_rate %>% 
              mutate(Position = str_remove_all(Position, "\\.0")),
            by = "Position") %>% 
  filter(Error == 0) %>% 
  group_by(Mom, Dad, Offspring) %>% 
  summarize(nb_SNPs = n(), .groups = "drop_last") %>%
  left_join(get_marker_types(data, "E4894", "E4733", "E5017")$Type_markers %>% 
              group_by(Mom, Dad) %>% 
              summarize(total_snps = n(), .groups = "drop_last"),
            by = c("Mom", "Dad")) %>% 
  mutate(Proportion = nb_SNPs / total_snps,
         Expected_prop_Mendel = case_when(
           Mom == 0 & Dad == 0 & Offspring == 0 ~ 1,
           Mom == 0 & Dad == 1 & Offspring != 2 ~ 0.5,
           Mom == 0 & Dad == 2 & Offspring == 1 ~ 1,
           Mom == 1 & Dad == 0 & Offspring != 2 ~ 0.5,
           Mom == 1 & Dad == 1 & Offspring != 1 ~ 0.25,
           Mom == 1 & Dad == 1 & Offspring == 1 ~ 0.5,
           Mom == 1 & Dad == 2 & Offspring != 0 ~ 0.5,
           Mom == 2 & Dad == 0 & Offspring == 1 ~ 1,
           Mom == 2 & Dad == 1 & Offspring != 0 ~ 0.5,
           Mom == 2 & Dad == 2 & Offspring == 2 ~ 1,
           TRUE ~ 0
         ),
         Expected_nb_SNPs = round(Expected_prop_Mendel * total_snps)) %>% 
  relocate(Mom, Dad, Offspring, nb_SNPs, Expected_nb_SNPs, total_snps, Proportion, Expected_prop_Mendel) %>% 
  filter(Mom == 2 & Dad == 0)


vcf <- read.vcfR("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/thinned_no_error_good_indivs_maf005_mindp20.vcf")
depth <- read.table("/shared/projects/sexisol/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Depth.tsv",
                    sep = "\t", header = TRUE)

# Forsmani
get_marker_types(data, "E4894", "E4733", "E5017")$Type_markers %>% 
  filter(Mom == 0 & Dad == 2) %>% 
  mutate(Position = str_remove_all(Position, "\\.1")) %>% 
  left_join(depth %>% 
              select(ID, E4894, E4733, E5017),
            by = join_by("Position" == "ID")) %>% 
  mutate(Depth_ref = str_split_fixed(E5017, ",", 2)[, 1] %>% as.numeric(),
         Depth_alt = (str_split_fixed(E5017, ",", 2)[, 2] %>% 
           str_split_fixed(., ":", 2))[, 1] %>% as.numeric(),
         Tot_depth = str_split_fixed(E5017, ":", 2)[, 2] %>% as.numeric(),
         Data = "Forsmani") %>% 
  select(-c(E4894, E4733, E5017)) %>% 
  rbind(get_marker_types(data, "E4907b", "E4748", "E5142")$Type_markers %>% 
          filter(Mom == 0 & Dad == 2) %>% 
          mutate(Position = str_remove_all(Position, "\\.1")) %>% 
          left_join(depth %>% 
                      select(ID, E4907b, E4748, E5142),
                    by = join_by("Position" == "ID")) %>% 
          mutate(Depth_ref = str_split_fixed(E5142, ",", 2)[, 1] %>% as.numeric(),
                 Depth_alt = (str_split_fixed(E5142, ",", 2)[, 2] %>% 
                                str_split_fixed(., ":", 2))[, 1] %>% as.numeric(),
                 Tot_depth = str_split_fixed(E5142, ":", 2)[, 2] %>% as.numeric(),
                 Data = "Prae") %>% 
          select(-c(E4907b, E4748, E5142))) %>% 
    rbind(get_marker_types(data, comb$Mom, comb$Dad, comb$Offspring)$Type_markers %>% 
            filter(Mom == 0 & Dad == 2) %>% 
            mutate(Position = str_remove_all(Position, "\\.1")) %>% 
            left_join(depth %>% 
                        select(ID, E5453, E4936, E4769),
                      by = join_by("Position" == "ID")) %>% 
            mutate(Depth_ref = str_split_fixed(E5453, ",", 2)[, 1] %>% as.numeric(),
                   Depth_alt = (str_split_fixed(E5453, ",", 2)[, 2] %>% 
                                  str_split_fixed(., ":", 2))[, 1] %>% as.numeric(),
                   Tot_depth = str_split_fixed(E5453, ":", 2)[, 2] %>% as.numeric(),
                   Data = "Hetero") %>% 
            select(-c(E5453, E4936, E4769))) %>% 
  pivot_longer(cols = contains("Depth"), names_to = "Allele", values_to = "depth") %>% 
    ggplot() +
    geom_boxplot(aes(x = factor(Offspring), y = depth, fill = Data)) +
    # facet_wrap(vars(Data))
  ggh4x::facet_grid2(Allele ~ Data) +
  my_theme

(depth %>% 
  select(ID, E5017, E5142, E5453) %>% 
  mutate(across(c(E5017, E5142, E5453), ~ str_split_fixed(., ",", 2)[, 1] %>% as.numeric(), .names = "{col}_ref"),
         across(c(E5017, E5142, E5453), ~ (str_split_fixed(., ",", 2)[, 2] %>% str_split_fixed(":", 2))[, 1] %>% as.numeric(), .names = "{col}_alt"),
         across(c(E5017, E5142, E5453), ~ str_split_fixed(., ":", 2)[, 2] %>% as.numeric(), .names = "{col}_tot")) %>% 
  pivot_longer(contains("_"), names_to = "Allele", values_to = "Depth") %>% 
  mutate(Data = case_when(str_split_fixed(Allele, "_", 2)[, 1] == "E5017" ~ "Forsmani",
                          str_split_fixed(Allele, "_", 2)[, 1] == "E5142" ~ "Praehirsuta",
                          TRUE ~ "Hetero"),
         Allele = str_split_fixed(Allele, "_", 2)[, 2] %>% str_to_title()) %>% 
  rename(loci = ID) %>%  
  prepare_data_for_manhattan() %>% 
  geom_manhattan(aes(y = Depth, facet1 = Allele, facet2 = Data))) +
  ggh4x::facet_grid2(Data ~ Allele)
    

# distrib <- 
depth %>% 
  select(ID, E5017, E5142, E5453) %>% 
  mutate(across(c(E5017, E5142, E5453), ~ str_split_fixed(., ",", 2)[, 1] %>% as.numeric(), .names = "{col}_ref"),
         across(c(E5017, E5142, E5453), ~ (str_split_fixed(., ",", 2)[, 2] %>% str_split_fixed(":", 2))[, 1] %>% as.numeric(), .names = "{col}_alt"),
         across(c(E5017, E5142, E5453), ~ str_split_fixed(., ":", 2)[, 2] %>% as.numeric(), .names = "{col}_tot")) %>% 
  pivot_longer(contains("_"), names_to = "Allele", values_to = "Depth") %>% 
  mutate(Data = case_when(str_split_fixed(Allele, "_", 2)[, 1] == "E5017" ~ "Forsmani",
                          str_split_fixed(Allele, "_", 2)[, 1] == "E5142" ~ "Praehirsuta",
                          TRUE ~ "Hetero"),
         Sample = str_split_fixed(Allele, "_", 2)[, 1],
         Allele = str_split_fixed(Allele, "_", 2)[, 2] %>% str_to_title())  %>% 
  filter(Depth < 100) %>% 
  # filter(log(Depth) > 6) %>% 
  ggplot(aes(x = Depth, fill = Allele)) +
  geom_histogram(bins = 25, alpha = 0.5, position = "dodge") +
  # geom_vline(aes(xintercept = log10(15)), colour = "red") +
  # geom_vline(aes(xintercept = log10(30)), colour = "red") +
  facet_wrap(vars(Data)) +
  my_theme

# boxplots <- 
depth %>% 
  select(ID, E5017, E5142, E5453) %>% 
  mutate(across(c(E5017, E5142, E5453), ~ str_split_fixed(., ",", 2)[, 1] %>% as.numeric(), .names = "{col}_ref"),
         across(c(E5017, E5142, E5453), ~ (str_split_fixed(., ",", 2)[, 2] %>% str_split_fixed(":", 2))[, 1] %>% as.numeric(), .names = "{col}_alt"),
         across(c(E5017, E5142, E5453), ~ str_split_fixed(., ":", 2)[, 2] %>% as.numeric(), .names = "{col}_tot")) %>% 
  pivot_longer(contains("_"), names_to = "Allele", values_to = "Depth") %>% 
  mutate(Data = case_when(str_split_fixed(Allele, "_", 2)[, 1] == "E5017" ~ "Forsmani",
                          str_split_fixed(Allele, "_", 2)[, 1] == "E5142" ~ "Praehirsuta",
                          TRUE ~ "Hetero"),
         Sample = str_split_fixed(Allele, "_", 2)[, 1],
         Allele = str_split_fixed(Allele, "_", 2)[, 2] %>% str_to_title()) %>% 
  filter(Depth < 100) %>% 
  ggplot(aes(x = Depth, y = Data, fill = Allele)) +
  geom_boxplot() +
  facet_wrap(vars(Allele)) +
  my_theme

depth %>% 
  select(ID, E5017, E5142, E5453) %>% 
  mutate(across(c(E5017, E5142, E5453), ~ str_split_fixed(., ",", 2)[, 1] %>% as.numeric(), .names = "{col}_ref"),
         across(c(E5017, E5142, E5453), ~ (str_split_fixed(., ",", 2)[, 2] %>% str_split_fixed(":", 2))[, 1] %>% as.numeric(), .names = "{col}_alt"),
         across(c(E5017, E5142, E5453), ~ str_split_fixed(., ":", 2)[, 2] %>% as.numeric(), .names = "{col}_tot")) %>% 
  pivot_longer(contains("_"), names_to = "Allele", values_to = "Depth") %>% 
  mutate(Data = case_when(str_split_fixed(Allele, "_", 2)[, 1] == "E5017" ~ "Forsmani",
                          str_split_fixed(Allele, "_", 2)[, 1] == "E5142" ~ "Praehirsuta",
                          TRUE ~ "Hetero"),
         Sample = str_split_fixed(Allele, "_", 2)[, 1],
         Allele = str_split_fixed(Allele, "_", 2)[, 2] %>% str_to_title()) %>% 
  group_by(Sample) %>% 
  summarize(quantile_25 = quantile(Depth, probs = 0.25, na.rm = TRUE),
            quantile_50 = quantile(Depth, probs = 0.50, na.rm = TRUE),
            quantile_75 = quantile(Depth, probs = 0.75, na.rm = TRUE),
            mean_depth = mean(Depth, na.rm = TRUE),
            Species = unique(Data))
  
ggpubr::ggarrange(distrib, boxplots, ncol = 1)
# Forsmani
get_marker_types(data, "E4894", "E4733", "E5017")$Summary_table %>% 
  filter(Mom == 0, Dad == 2)

# Prae
get_marker_types(data, "E4907b", "E4748", "E5142")$Summary_table %>% 
  filter(Mom == 0, Dad == 2)

# Heterospecific
get_marker_types(data, "E4935", comb$Dad, comb$Offspring)$Summary_table %>% 
  filter(Mom == 0, Dad == 2)

get_marker_types(data, comb$Mom, comb$Dad, comb$Offspring)$Type_markers %>% 
  filter(Mom == 2, Dad == 0, Offspring == 0) %>%
  add_markers_to_manhattan(Fst_prae_forsmani)

read.table("/shared/home/bpajot/Jaera/Softwares/Colony/small_data_with_mothers/small_data.PairwisePaternity",
           sep = ",", header = TRUE, comment.char = "|") %>% 
  mutate(prop = X.ExcLoci/X.PairGtype) %>% 
  ggplot(aes(x = prop)) +
  # geom_point()
  geom_histogram()
  
# 
# library(latticeExtra)
# 
# St <- Type_markers$Summary_table %>% 
#   ungroup() %>% 
#   mutate(across(c(Mom, Dad), ~ 2 - .),
#          across(c(Mom, Dad), ~ paste(., Offspring, sep = "_")),
#          across(c(Mom, Dad), ~ factor(., levels = c("0_0", "0_1", "0_2",
#                                                     "1_0", "1_1", "1_2",
#                                                     "2_0", "2_1", "2_2"))))
# 
# par(mfrow = c(1, 2))
# cloud(Proportion ~ Mom + Dad, St, panel.3d.cloud = panel.3dbars, col.facet='grey', 
#       xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
#       par.settings = list(axis.line = list(col = "transparent")))
# cloud(Expected_prop_Mendel ~ Mom + Dad, St, panel.3d.cloud = panel.3dbars, col.facet='grey', 
#       xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1), 
#       par.settings = list(axis.line = list(col = "transparent")))
# par(mfrow = c(1, 1))


error_rate <- read.table("/shared/home/bpajot/Jaera/finalresult/ddRAD_multiple_choice_exp/09_thin_vcf/Information_vcf/Error_rate_per_SNP.tsv",
           sep = " ", header = TRUE) %>% 
  filter(Error <= 0.1) %>% 
  ggplot(aes(x = Error)) + 
  geom_histogram()
