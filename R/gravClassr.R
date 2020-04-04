#' This function treats k means clusters as physical objects with mass and gravity,
#' calculating that gravity for each observation and flagging observations whose
#' distance and gravity classifications differ.
#'
#' This function requires three inputs:
#' 1. dataframe
#' 2. List of numeric variables for use in clustering algorithm
#' 3. Value of k
#'
#' This function produces one output, a dataframe with all relevant output.

gravClassr <- function(df, var, k) {

  # Load packages
  if (!require(tidyverse)) install.packages('tidyverse')
  if (!require(tidyselect)) install.packages('tidyselect')
  if (!require(glue)) install.packages('glue')
  library(tidyverse)
  library(tidyselect)
  library(glue)

  # Produce generic names for attributes for use in calculations
  newvar <- paste0("var", seq(1:length(var)))

  if(min(var %in% colnames(df)) == 0) {
    stop("One or more of the input columns does not exist in the dataframe")
  }
  else {

    # Subset df to only the clustering variables and rename
    x_df <-
      df %>%
      # variables from var
      select(all_of(var)) %>%
      # rename var to newvar
      rename(setNames(var, newvar))

    # Confirm columns are all numeric
    if(mean(sapply(x_df, is.numeric)) != 1) {
      stop("Parameter 'var' must contain only numeric columns")
    }
    else {

      # Handle missing cases
      if( sum(complete.cases(x_df)) != nrow(x_df) ) warning('Observations with missing rows will be removed')

      x_df_na_drop <-
        x_df %>%
        drop_na()

      # Call k-means algorithm
      x_result <- kmeans(x_df_na_drop, k)

      # Append cluster assignments
      x_df_na_drop$cluster <- x_result$cluster

      # Use row labels to name clusters
      x_centers <- as.data.frame(x_result$centers) %>%
        rownames_to_column('cluster')

      # Calculate mass and transpose
      x_mass <-
        x_df_na_drop %>%
        group_by(cluster) %>%
        summarize(mass = n()) %>%
        mutate(k = 1) %>%
        select(k, everything()) %>%
        gather(variable, value, -(k:cluster)) %>%
        unite(temp, variable, cluster) %>%
        spread(temp, value)

      # Gather centers and transpose
      x_centers_wide <-
        x_centers %>%
        mutate(k = 1) %>%
        select(k, everything()) %>%
        gather(variable, value, -(k:cluster)) %>%
        unite(temp, variable, cluster) %>%
        spread(temp, value)

      # Join mass and centers to df
      x_df_out <-
        x_df_na_drop %>%
        rownames_to_column('id') %>%
        mutate(id = as.numeric(id)) %>%
        mutate(k = 1) %>%
        inner_join(x_centers_wide,
                   by = 'k') %>%
        inner_join(x_mass,
                   by = 'k') %>%
        select(-k)

      # Calculate gravity
      for(i in 1:k) {
        varname <- paste0('grav_', i)
        dist_str <- glue_collapse(map_chr(.x = newvar, function(x) glue("({x} - {x}_{i})^2")), " + ")
        x_df_out <-
          mutate(x_df_out,
                 !!varname :=
                   eval(parse(text = paste0('mass_', i)))/eval(parse(text = dist_str)))
      }

      # Clean output
      x_df_out %<>%
        select(id, vars_select(names(.), cluster, starts_with('grav', ignore.case = TRUE))) %>%
        gather(col, val, vars_select(names(.), starts_with('grav', ignore.case = TRUE))) %>%
        group_by(id, cluster) %>%
        slice(which.max(val)) %>%
        mutate(grav_cluster = as.numeric(str_extract(col, "(\\d)+"))) %>%
        select(-col, -val)

      # Join to original data and create flag
      x_out <-
        df %>%
        rownames_to_column('id') %>%
        mutate(id = as.numeric(id)) %>%
        left_join(x_df_out,
                   by = 'id') %>%
        mutate(flag = if_else(cluster == grav_cluster, 0, 1)) %>%
        bind_rows(x_centers %>%
                    rename(setNames(newvar, var)) %>%
                    rownames_to_column('id') %>%
                    mutate(id = nrow(df) + as.numeric(id),
                           flag = 2,
                           cluster = as.numeric(cluster)))

        return(x_out)

    }
  }
}
