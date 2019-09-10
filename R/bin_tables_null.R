bin_tables_null <- function(M_folder, M_format = "shp", vars_folder,
                            var_format = "GTiff", round = FALSE, round_names,
                            multiplication_factor = 1, bin_size = 10,
                            output_folder = "Null_species_E_space_bins"){
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

  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_folder)

  bin_tabs <- lapply(1:dim(variables)[3], function(i) {
    cat("\nPreparing range values from environmental layers:\n")
    M_range <- lapply(1:length(mlist), function(j) {
      if (M_format %in% c("shp", "gpkg")) {
        if (M_format == "shp") {
          M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j], verbose = FALSE)
        } else {
          M <- rgdal::readOGR(paste0(M_folder, "/", mlist[j]), spnames[j], verbose = FALSE)
        }
      } else {
        M <- raster::raster(mlist[j])
      }

      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      range <- range(na.omit(mvar[]))
      cat("\t", j, "of", length(mlist), "species finished\n")

      return(range)
    })

    # overall range
    M_range <- do.call(rbind, M_range)
    overall_range <- range(c(M_range))

    if (overall_range[2] > 999) {
      overall_range <- round(overall_range / 10)
      M_range <- round(M_range / 10)
    }
    if (overall_range[2] > 9999) {
      overall_range <- round(overall_range / 100)
      M_range <- round(M_range / 100)
    }

    # modification of range
    o_minimum <- overall_range[1]
    o_minimumc <- ifelse(o_minimum == 0, 0, floor(o_minimum / bin_size) * bin_size) - bin_size

    o_maximum <- overall_range[2]
    o_maximumc <- ifelse(o_maximum == 0, 0, ceiling(o_maximum / bin_size) * bin_size) + bin_size

    overall_range <- c(o_minimumc, o_maximumc)

    # bin tables
    cat("\nPreparing bin tables using ranges:\n")

    bin_table <- bin_env_null(overall_range, M_range, bin_size)
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
