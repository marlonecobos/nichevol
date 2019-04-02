bins_table <- function(M_folder, M_format = "shp", occ_folder, longitude_col, latitude_col,
                       vars_folder, vars_format = "GTiff", round = FALSE, round_names,
                       multiplication_factor = 1, percentage_out = 5, bin_size = 10,
                       output_folder = "Species_E_space_bins"){
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

  for (i in 1:dim(variables)[3]) {
    # data
    M_range <- list()
    sp_range <- list()

    cat("\nPreparing range values from environmental layers and species data:\n")

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
      df_layer <- abs(mval - medians)
      names(df_layer) <- mval

      ## limit
      limit <- floor((100 - percentage_out) * length(df_layer) / 100)
      df_layer <- df_layer[1:limit]

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

    o_minimum <- overall_range[1]
    o_minimumc <- ifelse(o_minimum == 0, 0, floor(o_minimum / 10) * 10) - 10

    o_maximum <- overall_range[2]
    o_maximumc <- ifelse(o_maximum == 0, 0, ceiling(o_maximum / 10) * 10) + 10
    overall_range <- c(o_minimumc, o_maximumc)

    # bin tables
    cat("\nPreparing bin tables using ranges:\n")

    bin_table <- bins_env(overall_range, M_range, sp_range, bin_size)
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(gsub("_", " ", spnames), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    # write table
    write.csv(bin_table, paste0(output_folder, "/", names(variables)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "OF", dim(variables)[3], "VARIABLES FINISHED\n")
  }
}


bins_env <- function(overall_range, M_range, sp_range, bin_size){
  # sequences
  sequence_vals <- seq(overall_range[1], overall_range[2], bin_size)

  # numeric ranges
  ranges <- sapply(1:(length(sequence_vals) - 1), function(x){
    if (x == 1) {
      c(sequence_vals[x], sequence_vals[x + 1])
    } else {
      c(sequence_vals[x] + 1 , sequence_vals[x + 1])
    }
  })

  # character bins
  bins <- sapply(1:(length(sequence_vals) - 1), function(x){
    if (x == 1) {
      paste(sequence_vals[x], sequence_vals[x + 1], sep = " to ")
    } else {
      paste(sequence_vals[x] + 1 , sequence_vals[x + 1], sep = " to ")
    }
  })

  bin_tab <- list()

  for (i in 1:dim(M_range)[1]) {
    M_test <- seq(M_range[i, 1], M_range[i, 2], 1)
    sp_test <- seq(sp_range[i, 1], sp_range[i, 2], 1)

    invar_M <- vector()
    invar_sp <- vector()

    for (j in 1:dim(ranges)[2]) {
      if (sum(M_test %in% seq(ranges[1, j], ranges[2, j], 1)) > 0) {
        invar_M[j] <- 1
      } else {
        invar_M[j] <- 0
      }

      if (sum(sp_test %in% seq(ranges[1, j], ranges[2, j], 1)) > 0) {
        invar_sp[j] <- 100
      } else {
        invar_sp[j] <- 0
      }
    }

    invar_sum <- invar_M + invar_sp

    bin_tab[[i]] <- vector()
    for (j in 1:length(invar_sum)) {
      if(invar_sum[j] == 0) bin_tab[[i]][j] <- "?"
      if(invar_sum[j] == 1) bin_tab[[i]][j] <- "0"
      if(invar_sum[j] >= 100) bin_tab[[i]][j] <- "1"
    }

    cat("\t", i, "of", dim(M_range)[1], "species finished\n")
  }

  bin_tab <- do.call(rbind, bin_tab)
  colnames(bin_tab) <- bins

  return(bin_tab)
}
