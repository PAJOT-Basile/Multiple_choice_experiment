#' This script contains transversal useful functions that I use in my code and
#' that I do not want to re-define in each script. In addition, here, the functions
#' are commented and a docstring is added, hopefully making the functions a little
#' bit more comprehensible (sorry about some of the functions that are a horror).


#################################
###### Random functions #########
#################################
"%not_in%" <- function(x, y){
  #' Not in
  #' 
  #' This function is the opposite of the "%in%" function. It checks if the 
  #' values contained in input x are NOT contained in y.
  #' 
  #'@param x (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the y vector
  #'@param y (vector).
  #'    This input contains any type of object, as long as they are the same type
  #'    as the ones in the x vector
  #'    
  #'@returns (bool)
  #'    This function returns TRUE for each value of x that is NOT in y and FALSE
  #'    for every value of x that is in y.
  #' @export
  
  return(!(x %in% y))
}

str_split_last <- function(string, pattern){
  #' Split string on last pattern
  #' 
  #' @description
    #' This function works like the `str_split` functions, that allows to cut
    #' a string on the last occurrence of a pattern in the string
  #' 
  #' @param string String to cut
  #' @param pattern Pattern on which to cut string
  #' 
  #' @return A string that contains all the provided string except what comes after
  #' the last occurrence of the pattern in the string.
  nb_substrings <- str_count(string, pattern)
  nb_substrings <- ifelse(pattern != "", nb_substrings + 1, nb_substrings)
  return(str_split_fixed(string, pattern, nb_substrings)[, nb_substrings])
}

#################################
###### Summary functions ########
#################################
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


#################################
###### Popgen functions #########
#################################
compute_Hobs <- function(genetic_data){
  #' Compute observed heterozygosity
  #' 
  #' @param genetic_data Genind object. This is the genind object to use on which
  #' to compute the heterozygosity
  #' 
  #' @returns A named vector. The values of the vector are the observed heterozygosity
  #' for each position in the genind object and the names are the names of the
  #' position in the genind object
  NA_count <- colSums(is.na(genetic_data@tab))
  N_het <- colSums(genetic_data@tab == 1, na.rm = TRUE)
  N_geno <- length(genetic_data@tab[, 1]) - NA_count
  Hobs <- N_het / N_geno
  return(Hobs)
}
#################################
###### Plotting functions #######
#################################
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
  
  geom_box_background <- function(To_plot_manhat, colName, chromosome_centers){
    #' Make polygon delimitations for geom_polygon
    #' 
    #' This function creates a data frame that contains the polygon delimitations
    #' if it is required to separate chromosomes by using boxes in the visual
    #' representations
    #' 
    #'@param To_plot_manhat (data.frame).
    #'    This data frame contains the data on the plotting of the values in a 
    #'    manhattan plot. 
    #'@param colName (string or name).
    #'    This argument is the name of the column in "To_plot_manhat" that we want
    #'    to plot.
    #'@param chromosomes_centers (data.frame).
    #'    This data frame contains information on the chromosomes: their name and
    #'    the position of their center is the minimum required for this function to
    #'    work.
    #'    
    #'@returns (data.frame)
    #'    This data frame contains the information on delimitations of the polygons
    #'    (position of the angles on the x and y axis). 
    #' @export
    
    # First, we transform the colName to use it with tidyverse
    colName <- as.character(colName)
    colName <- as.name(substitute(colName))
    
    # Then, we make y delimitations (max and min) for the boxes. The height of the
    # boxes has to be bigger than the maximal value and minimal value of the
    # considered column so we add/remove 20% to its value
    max_col <- (To_plot_manhat %>% 
                  select(as.symbol(colName)) %>% max(na.rm = TRUE)) * 1.2
    min_col <- To_plot_manhat %>% 
      select(as.symbol(colName)) %>% min(na.rm = TRUE)
    # Here, we just consider the case where the value of the minimum is positive
    # or negative
    min_col <- min_col * ifelse(min_col >= 0, 0.8, 1.2)
    
    
    # Then, we add this information as well as the x delimitations of the
    # chromosomes to the data frame to return 
    boxes <- To_plot_manhat  %>% 
      # First, we merge the data with the data frame containing the center of the
      # chromosomes to acess their center position
      right_join(chromosome_centers, by="Chromosome") %>% 
      # We group the data frame by the chromosome to treat each chromosome
      # separately
      group_by(Chromosome) %>% 
      # Then, we calculate the max and min positions on each chromosome
      summarise(Min = bp_cum %>% min,
                Max = bp_cum %>% max) %>% 
      # As the geom_polygon function requires to have all the x values in the same
      # column, we pivot the data frame to regroup the max and min value per 
      # chromosome
      pivot_longer(cols = c(Max, Min), values_to = "x_value", names_to = "operation") %>% 
      # The column "operation" contains the names max an min, so we get rid of it
      # because it does not bring any new information
      select(-operation) %>% 
      # As the geom_polygon function requires to have the positions of all the
      # polygon angles on both x and y axis, this means that we need to have the 
      # x/y positions of the first angle, followed by the one next to it, and so
      # on until you run out of angles. So, as we are building a rectangle, we need
      # 4 angles, so 4 positions along the x axis and 4 positions along the y axis.
      # Therefore, we duplicate the lines we already created for each chromosome
      rbind(., .) %>% 
      # We arrange the data frame by Chromosome and x value
      arrange(Chromosome, x_value) %>% 
      # And finally, we add the y values that need to have this order to be able
      # to build the polygons.
      cbind("y_value" = rep(c(min(c(0, min_col)),
                              max(c(1.02, max_col)),
                              max(c(1.02, max_col)),
                              min(c(0, min_col))),
                            nrow(.)/4))
    
    # Return the created variables
    return(boxes)
  }
  
  is.continuous <- function(x){
    #' Checks for continuousity
    #' 
    #' This function checks if a variable is discrete or continuous.
    #' 
    #' @param x (numeric vector).
    #'      A vector containing values to see if it is continuous or discrete
    #' @returns (bool).
    #'      This function returns TRUE if the variable is continuous (more than 10
    #'      levels) or FALSE if the variable is discrete (less than 10 levels)
    #' @section Warning:
    #' Here a variable is considered discrete if it contains less than 10
    #' distinct levels. This approximation is done to simplify the distinction 
    #' between continuous or discrete variables, but it is not a real classification
    #' @export
    
    return(length(unique(x)) >= 10)
  }
  
  # First, we added some error checking to keep the function from running for
  # nothing
  if ("Position" %not_in% names(df)){
    stop("This function needs a column called 'Position' that contains the position of each SNP on the genome")
  }
  if (is.numeric(df$Position[1])){
    if ("Chromosome" %not_in% names(df)){
      stop("If the positions are already usable, please name a column 'Chromosome' to be used as a chromosome reference.")
    }
  }
  if ("y" %not_in% names(mapping)){
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
               factor(levels = df %>% select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
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
                        select(Chromosome) %>% 
                        unique %>% 
                        arrange(as.numeric(gsub("\\D*(\\d+).*", "\\1", Chromosome))) %>% 
                        as.vector %>% unname %>% unlist)) %>% 
      # And we select the columns of interest
      select(Position, Chromosome, matches(to_select)) %>% 
      # Finally, we drop the missing values
      drop_na()
  }
  
  # Then, we make a cumulative data frame with the positions of each end and beginning of chromosome
  data_cum <- To_plot %>% 
    group_by(Chromosome) %>% 
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
    if (param %not_in% names(mapping)){
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



