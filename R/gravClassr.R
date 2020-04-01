# This function treats k means clusters as physical objects with mass and gravity,
# calculating that gravity for each observation and flagging observations whose
# distance and gravity classifications differ.
#
# This function requires three inputs:
# 1. dataframe
# 2. List of numeric variables for use in clustering algorithm
# 3. Value of k
#
# This function produces one output, a dataframe with all relevant output.

gravClassr <- function(df, var, k) {

  # Produce generic names for attributes for use in calculations
  newvar <- paste0("var", seq(1:length(var)))

  # Subset df to only the clustering variables and rename
  x_df <-
    df %>%
    # variables from var
    select(all_of(var)) %>%
    # rename var to newvar
    rename(setNames(var, newvar))

  # Call k-means algorithm
  x_result <- kmeans(x_df, k)

  # Append cluster assignments
  x_df$cluster <- x_result$cluster

  # Use row labels to name clusters
  x_centers <- as.data.frame(x_result$centers) %>%
    rownames_to_column('cluster')

  # Calculate mass and transpose
  x_mass <-
    x_df %>%
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
    x_df %>%
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
    select(vars_select(names(.), cluster, starts_with('grav', ignore.case = TRUE))) %>%
    rownames_to_column('id') %>%
    mutate(id = as.numeric(id)) %>%
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
    inner_join(x_df_out,
               by = 'id') %>%
    mutate(flag = if_else(cluster == grav_cluster, 0, 1))

  return(x_out)
}
