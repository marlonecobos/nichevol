#' Bin tables of environemntal conditions in M and occurrences
#'
#' @description bin_tables0 helps in creating csv files with bin tables
#' of environmental conditions in M and species occurrence records. All of this
#' starting from raw data, for various species, and using multiple variables.
#'
#' @param M_folder (character) name of the folder containing files representing
#' the accessible area (M) for all species to be analyzed.
#' @param M_format format of files representing the accessible area (M) for the
#' species. Names must match with the ones for occurrence files in
#' \code{occ_folder}. Options are: "shp", "gpkg", or any of the options in
#' \code{\link[raster]{writeFormats}}.
#' @param occ_folder (character) name of the folder containing csv files of
#' occurrence data for all species.
#' @param longitude (character) name of the column containing values of longitude
#' of occurrences.
#' @param latitude (character) name of the column containing values of latitude
#' of occurrences.
#' @param var_folder (character) name of the folder conatining layers to
#' represent environmental variables.
#' @param var_format format of layers to represent environmental variables. See
#' options in \code{\link[raster]{writeFormats}}.
#' @param round (logical) whether or not to round the values of one or more
#' variables after multiplying them times \code{multiplication_factor}.
#' Default = FALSE.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. Default = 5.
#' @param bin_size (numeric) size of bins. Interval of values to be considered
#' when creating bin tables. Default = 10.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Species_E_bins".
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all resultant csv
#' files for all variables, with bins for all species.

bin_tables0 <- function(M_folder, M_format, occ_folder, longitude,
                        latitude, var_folder, var_format,
                        round = FALSE, round_names, multiplication_factor = 1,
                        percentage_out = 5, bin_size = 10,
                        output_directory = "Species_E_bins"){
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument M_folder is missing.")}
  if (missing(M_format)) {stop("Argument M_format is missing.")}
  if (missing(occ_folder)) {stop("Argument occ_folder is missing.")}
  if (missing(longitude)) {stop("Argument longitude is missing.")}
  if (missing(latitude)) {stop("Argument latitude is missing.")}
  if (missing(var_folder)) {stop("Argument var_folder is missing.")}
  if (missing(var_format)) {stop("Argument var_format is missing.")}

  # formats and data to start
  cat("\nPreparing data, please wait...\n\n")
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

  occlist <- list.files(path = occ_folder, pattern = ".csv$", full.names = TRUE)
  variables <- raster::stack(list.files(path = var_folder, pattern = v_patt,
                                        full.names = TRUE))

  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_directory)
  cat("Preparing range values and bin tables from environmental layers and species data:\n")
  bin_tabs <- lapply(1:dim(variables)[3], function(i) {
    # data
    M_range <- list()
    sp_range <- list()

    cat("\n   Preparing range values:\n")

    for (j in 1:length(occlist)) {
      ## M
      if (M_format %in% c("shp", "gpkg")) {
        if (M_format == "shp") {
          M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j], verbose = FALSE)
        } else {
          M <- rgdal::readOGR(paste0(M_folder, "/", mlist[j]), spnames[j],
                              verbose = FALSE)
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
      occval <- na.omit(raster::extract(mvar, occ[, c(longitude,
                                                      latitude)]))
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
    o_minimumc <- ifelse(o_minimum == 0, 0,
                         floor(o_minimum / bin_size) * bin_size) - bin_size

    o_maximum <- overall_range[2]
    o_maximumc <- ifelse(o_maximum == 0, 0,
                         ceiling(o_maximum / bin_size) * bin_size) + bin_size

    overall_range <- c(o_minimumc, o_maximumc)

    # bin tables
    cat("   Preparing bin tables using ranges:\n")

    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(gsub("_", " ", spnames), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    # write table
    write.csv(bin_table,
              paste0(output_directory, "/", names(variables)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "of", dim(variables)[3], "variables processed\n")
    return(bin_table)
  })

  names(bin_tabs) <- names(variables)
  return(bin_tabs)
}
