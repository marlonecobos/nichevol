histograms_env <- function(M_folder, M_format = "shp", occ_folder, longitude_col, latitude_col,
                           vars_folder, vars_format = "GTiff", CL_lines = c(95, 99), col,
                           round = FALSE, round_names, multiplication_factor = 1,
                           output_folder = "Histogram_ranges_check", save_ranges = TRUE){
   # checking for potential errors
  if (missing(col)) {
    col <- sort(gray.colors(length(CL_lines) + 1), decreasing = TRUE)[1:length(CL_lines)]
  }

  # formats
  if (M_format == "shp") {
    M_patt <- ".shp$"
    subs <- ".shp"
  }
  if (M_format == "ascii") {
    M_patt <- ".asc$"
    subs <- ".asc"
  }
  if (M_format == "GTiff") {
    M_patt <- ".tif$"
    subs <- ".tif"
  }
  if (M_format == "EHdr") {
    M_patt <- ".bil$"
    subs <- ".tif"
  }
  if (vars_format == "ascii") {v_patt <- ".asc$"}
  if (vars_format == "GTiff") {v_patt <- ".tif$"}
  if (vars_format == "EHdr") {v_patt <- ".bil$"}

  # needed data to start
  cat("\nPreparing data, please wait...\n\n")

  mlist <- gsub(subs, "",list.files(path = M_folder, pattern = paste0(".", M_format, "$")))
  occlist <- list.files(path = occ_folder, pattern = ".csv$", full.names = TRUE)
  spnames <- gsub(".csv", "", list.files(path = occ_folder, pattern = ".csv$"))
  variables <- raster::stack(list.files(path = vars_folder, pattern = v_patt, full.names = TRUE))

  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_folder)

  ranges <- list()

  for (i in 1:dim(variables)[3]) {
    # data
    df_layer <- list()
    occ_dfs <- list()

    M_ranges <- list()
    M_limits <- list()
    sp_ranges <- list()
    y_values <- list()

    for (j in 1:length(occlist)) {
      ## M shapefiles
      M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j], verbose = FALSE)

      ## occurrences
      occ <- read.csv(occlist[j])

      # processing
      ## get values of variables in M
      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      mval <- na.omit(raster::values(mvar))

      ## distance of each absolute value to median value
      medians <- median(mval)
      df_layer[[j]] <- abs(mval - medians)
      names(df_layer[[j]]) <- mval

      occval <- na.omit(raster::extract(mvar, occ[, c(longitude_col, latitude_col)]))
      occ_dfs[[j]] <- abs(occval - medians)
      names(occ_dfs[[j]]) <- occval

      y_values[[j]] <- c(max(table(mval)), max(table(df_layer[[j]])))

      ## ranges of real values
      M_limit <- list()
      for (k in 1:length(CL_lines)) {
        limit <- floor((CL_lines[k]) * length(df_layer[[j]]) / 100)
        df_layera <- sort(df_layer[[j]])[1:limit]

        M_limit[[k]] <- range(as.numeric(names(df_layera)))
      }

      M_ranges[[j]] <- range(mval)
      M_limits[[j]] <- do.call(c, M_limit)
      sp_ranges[[j]] <- range(occval)

      cat("\t", j, "of", length(occlist), "species finished\n")
    }

    ## ranges final values for variables
    limits <- do.call(rbind, M_limits)
    ranges[[i]] <- data.frame(gsub("_", " ", spnames), do.call(rbind, sp_ranges),
                              do.call(rbind, M_ranges), limits)
    colnames(ranges[[i]]) <- c("Species", "Species_lower", "Species_upper", "M_lower", "M_upper",
                            paste0("M_", rep(CL_lines, each = 2), c("_lowerCL", "_upperCL")))

    if (save_ranges == TRUE) {
      write.csv(ranges[[i]], file = paste0(output_folder, "/Ranges_", names(variables)[i], ".csv"),
                row.names = FALSE)
    }

    ## frecuency of each value in M
    pdf_histograms(env_data = df_layer, occ_data = occ_dfs, y_values = y_values, sp_names = spnames,
                   variable_name = names(variables)[i], CL_lines = CL_lines, limits = limits,
                   col = col, output_folder = output_folder)

    cat(i, "of", dim(variables)[3], "variables processed\n")
  }

  if (save_ranges == TRUE) {
    cat("\nsave_ranges = TRUE:\ncsv files with the environmental ranges were saved in", output_folder, "\n")
  }

  names(ranges) <- names(variables)

  return(ranges)
}
