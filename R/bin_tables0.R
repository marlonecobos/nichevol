bin_tables0 <- function(M_folder, M_format = "shp", occ_folder, longitude_col, latitude_col,
                        vars_folder, vars_format = "GTiff", round = FALSE, round_names,
                        multiplication_factor = 1, percentage_out = 5, bin_size = 10,
                        output_folder = "Species_E_space_bins"){
  # formats and data to start
  if (M_format %in% c("shp", "gpkg")) {
    if (M_format == "shp") {
      M_patt <- ".shp$"
      subs <- ".shp"
      mlist <- gsub(subs, "", list.files(path = M_folder, pattern = M_patt))
      spnames <- mlist
    } else {
      M_patt <- ".gpkg$"
      subs <- ".gpkg"
      mlist <- list.files(path = M_folder, pattern = M_patt)
      spnames <- gsub(subs, "", mlist)
    }
  } else {
    M_patt <- paste0(rformat_type(var_format), "$")
    subs <- rformat_type(var_format)
    mlist <- list.files(path = M_folder, pattern = M_patt, full.names = TRUE)
    spnames <- gsub(subs, "", list.files(path = M_folder, pattern = M_patt))
  }
  v_patt <- paste0(rformat_type(var_format), "$")

  variables <- raster::stack(list.files(path = vars_folder, pattern = v_patt, full.names = TRUE))
  occlist <- list.files(path = occ_folder, pattern = ".csv$", full.names = TRUE)

  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_folder)

  bin_tabs <- lapply(1:dim(variables)[3], function(i) {
    # data
    M_range <- list()
    sp_range <- list()

    cat("\nPreparing range values from environmental layers and species data:\n")

    for (j in 1:length(occlist)) {
      ## M
      if (M_format %in% c("shp", "gpkg")) {
        if (M_format == "shp") {
          M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j], verbose = FALSE)
        } else {
          M <- rgdal::readOGR(paste0(M_folder, "/", mlist[j]), spnames[j], verbose = FALSE)
        }
      } else {
        M <- raster::raster(mlist[j])
      }

      ## occurrences
      occ <- read.csv(occlist[j])

      # processing
      ## get values of variables in M
      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      mval <- na.omit(mvar[])

      ## distance of each absolute value to median value
      medians <- median(mval)
      df_layer <- abs(mval - medians)
      names(df_layer) <- mval

      ## limit
      limit <- floor((100 - percentage_out) * length(df_layer) / 100)
      df_layer <- sort(df_layer)[1:limit]

      M_range[[j]] <- range(as.numeric(names(df_layer)))

      ## occurrences
      occval <- na.omit(raster::extract(mvar, occ[, c(longitude_col, latitude_col)]))
      sp_range[[j]] <- range(occval)

      cat("\t", j, "of", length(occlist), "species finished\n")
    }

    # overall range
    M_range <- do.call(rbind, M_range)
    sp_range <- do.call(rbind, sp_range)
    overall_range <- range(c(c(M_range), c(sp_range)))

    if (overall_range[2] > 999) {
      overall_range <- round(overall_range / 10)
      M_range <- round(M_range / 10)
      sp_range <- round(sp_range / 10)
    }
    if (overall_range[2] > 9999) {
      overall_range <- round(overall_range / 100)
      M_range <- round(M_range / 100)
      sp_range <- round(sp_range / 100)
    }

    # modification of range
    o_minimum <- overall_range[1]
    o_minimumc <- ifelse(o_minimum == 0, 0, floor(o_minimum / bin_size) * bin_size) - bin_size

    o_maximum <- overall_range[2]
    o_maximumc <- ifelse(o_maximum == 0, 0, ceiling(o_maximum / bin_size) * bin_size) + bin_size

    overall_range <- c(o_minimumc, o_maximumc)

    # bin tables
    cat("\nPreparing bin tables using ranges:\n")

    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(gsub("_", " ", spnames), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    # write table
    write.csv(bin_table, paste0(output_folder, "/", names(variables)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "of", dim(variables)[3], "variables processed\n")
    return(bin_table)
  })

  names(bin_tabs) <- names(variables)
  return(bin_tabs)
}
