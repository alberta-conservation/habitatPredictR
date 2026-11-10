# Function for converting state class, age, harvest rasters into dataframes for predicting species distributions with
# the ABMI species models via the ABMIexploreR package.

#' Title
#'
#' @param reference A True/FALSE statement
#' @param states A raster of state classes
#' @param ages A raster of ages
#' @param grid A polygon file of the
#' @param harvest A raster of cumulative harvest
#' @param seismic A raster of seismic line occurrence
#'
#' @returns A data frame
#' @export
#'
#' @examples
make_habitat_grid <- function(reference = TRUE,
                        states = NULL,
                        ages = NULL,
                        grid = NULL,
                        harvest = NULL,
                        seismic = NULL){

  # Which variable set to use (vectors embedded)
  if(reference){
    vars <- ref_variables
  }else{
    vars <- current_variables
  }

  ## Convert treed state classes into state class and age, including clear cuts if using harvest
  # Treed classes that included age age classes in the model (embedded)
  treed_ages <- c("Decid", "Mixedwood", "Pine", "Spruce", "TreedBog", "TreedFen", "TreedSwamp")
  # List of treed age classes (embedded)
  treed_ageClass <- c("R", 1:8)
  # Treed classes that can be designated as clear cuts in the model (embedded)
  cc_ages <- c("Decid", "Mixedwood", "Pine", "Spruce")
  # Age classes used for clear cuts (embedded)
  cc_ageClass <- c("R", 1:4)

  # Convert continuous ages into age classes
  age <- values(ages)
  av <- ifelse(age > 80, 80, age) # the maximum age class in abmi habitat models is 80
  age_class = as.character(cut(av,
                               breaks = seq(0, 80, 10),
                               labels = as.character(seq(1, 8, 1) - 1),
                               right = TRUE))
  age_class = if_else(av[, 1] == 80, "8", as.character(age_class))
  age_class = data.frame(age_class = if_else(age_class == 0, "R", age_class))

  ageClass <- ages
  values(ageClass) <- age_class$age_class
  names(ageClass) <- "hab_age"

  ##  Convert the state class map to model habitat classes
  # Create the lookup table
  habclass_lookup <- state_class_hf_lookup |> distinct(Description, .keep_all = TRUE)
  # Convert the state class raster to abmi habitat classes
  habClass <- terra::subst(x = states, from = state_class_hf_lookup$Id, to = state_class_hf_lookup$hab.id)
  # Convert habClass from numeric to factor
  levels(habClass) <- data.frame(value = habclass_lookup$hab.id, hab_class = habclass_lookup$Description)

  # Create the list of habitat proportions for each grid cell
  # Create the empty dataset and fill from each list item

  hP_list <- lapply(1:nrow(grid), function(i){
    st <- crop(habClass, vect(grid[i, ]), mask = TRUE) # cell-level habitat class
    ag <- crop(ageClass, vect(grid[i, ]), mask = TRUE) # cell-level age data
    if(!is.null(harvest)){
      hv <- crop(harvest, vect(grid[i, ]), mask = TRUE)
    }

    if(!is.null(seismic)){
      sm <- data.frame(sm = values(crop(seismic, vect(grid[i, ]), mask = TRUE))[, 1])
    }

    # Get the state and age data, and create the habitat-age classes for treed classes with age
    if(is.null(harvest)){
      grStateAge <- data.frame(state = values(st, dataframe = TRUE),
                               age = values(ag, dataframe = TRUE)) |>
        mutate(habitat = if_else(hab_class %in% treed_ages, paste0(hab_class, hab_age), hab_class),
               coef = landscape_lookup$COEF[match(habitat, landscape_lookup$ID)])
    }else{
      grStateAge <- data.frame(state = values(st, dataframe = TRUE),
                               age = values(ag),
                               values(hv)) |>
        rename(harv = names(hv)) |>
        mutate(habitat = if_else(hab_class %in% treed_ages, paste0(hab_class, hab_age), hab_class)) |>
        mutate(habitat = if_else(hab_class %in% cc_ages & hab_age %in% c(1:4, 9) & harv == 2, paste0("CC", hab_class, hab_age), habitat),
               coef = landscape_lookup$COEF[match(habitat, landscape_lookup$ID)])
    }

    if(!is.null(seismic)){
      grStateAge <- grStateAge |>
        mutate(seismic = sm$sm) |>
        mutate(seismic = ifelse(seismic == 1, 1 - median(grid_conv$prop), 1 - median(grid_low$prop))) |>
        mutate(seismic = ifelse(is.na(seismic), 1, seismic))
    }

    # Convert amounts into proportions
    if(!is.null(seismic)){
      hp <- drop_na(grStateAge) |>
        group_by(coef) |>
        summarise(prop = sum(seismic)/nrow(drop_na(grStateAge)))
      hp1 <- data.frame(t(hp |> select(prop))) |>
        set_names(hp$coef) |>
        mutate(EnSeismic = 1 - sum(hp$prop))
    }else{
      hp1 <- as.data.frame(as.list(table(grStateAge$coef)/nrow(grStateAge)))
    }

    d <- setNames(data.frame(matrix(ncol = length(vars), nrow = 0)), vars)
    d <- add_row(d, hp1)
    d[is.na(d)] <- 0
    rownames(d) <- grid$fid[i]

    return(d)
  })

  # Create the empty dataset and fill from each list item
  habProps <- bind_rows(hP_list)

  return(habProps)

}




