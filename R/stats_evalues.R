#' Statistics of environmental conditions in M and occurrences
#'
#' @description stats_evalues helps in creating csv files with statistics
#' of environmental conditions in M and species occurrence records, starting
#' from raw data for multiple species and using multiple variables.
#'
#' @param stats (character) name or vector of names of functions to be applied
#' to get basic statistics of environmental values.
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
#' @param var_folder (character) name of the folder containing layers to
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
#' @param output_directory (character) name of the folder in which results will
#' be written. Default = "Species_E_stats".
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all resultant csv
#' files for all variables, with statistics for all species. Results will also be
#' returned as a list.

stats_evalues <- function(stats = c("median", "range"), M_folder, M_format,
                          occ_folder, longitude, latitude, var_folder,
                          var_format, round = FALSE, round_names,
                          multiplication_factor = 1, percentage_out = 0,
                          output_directory = "Species_E_stats") {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}

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
  variables <- raster::stack(list.files(path = var_folder,
                                        pattern = v_patt, full.names = TRUE))

  var_names <- names(variables)
  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    noround <- var_names[!var_names %in% round_names]
    variables <- raster::stack(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_directory)
  cat("Preparing statistics from environmental layers and species data:\n")
  n_vars <- raster::nlayers(variables)
  var_stats <- lapply(1:n_vars, function(i) {
    sp_stats <- lapply(1:length(occlist), function(j) {
      ## data
      if (M_format %in% c("shp", "gpkg")) {
        if (M_format == "shp") {
          M <- rgdal::readOGR(dsn = M_folder, layer = mlist[j],
                              verbose = FALSE)
        } else {
          M <- rgdal::readOGR(paste0(M_folder, "/",
                                     mlist[j]), spnames[j], verbose = FALSE)
        }
      } else {
        M <- raster::raster(mlist[j])
      }
      occ <- read.csv(occlist[j])

      ## preparing e values
      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      mval <- na.omit(mvar[])
      if (percentage_out > 0) {
        medians <- median(mval)
        df_layer <- abs(mval - medians)
        names(df_layer) <- mval
        limit <- floor((100 - percentage_out) * length(df_layer)/100)
        df_layer <- sort(df_layer)[1:limit]
        mval <- as.numeric(names(df_layer))
      }
      occval <- na.omit(raster::extract(mvar, occ[, c(longitude, latitude)]))

      ## obtaining statistics
      if (length(stats) > 1) {
        m_stats <- lapply(1:length(stats), function(k) {
          eval(parse(text = paste0(stats[k], "(mval)")))
        })
        o_stats <- lapply(1:length(stats), function(k) {
          eval(parse(text = paste0(stats[k], "(occval)")))
        })

      } else {
        m_stats <- eval(parse(text = paste0(stats, "(mval)")))
        o_stats <- eval(parse(text = paste0(stats, "(occval)")))
      }
      names(m_stats) <- stats
      names(o_stats) <- stats

      cat("\t", j, "of", length(occlist), "species finished\n")
      return(list(M = unlist(m_stats), Occurrences = unlist(o_stats)))
    })

    # preparing tables with results
    m_table <- data.frame(Species = mlist,
                          do.call(rbind, lapply(sp_stats, function(x) {x[[1]]})))
    o_table <- data.frame(Species = mlist,
                          do.call(rbind, lapply(sp_stats, function(x) {x[[2]]})))

    # write table
    write.csv(m_table,
              paste0(output_directory, "/", var_names[i], "_M_stats.csv"),
              row.names = FALSE)
    write.csv(o_table,
              paste0(output_directory, "/", var_names[i], "_Occurrence_stats.csv"),
              row.names = FALSE)

    cat(i, "of", n_vars, "variables processed\n")

    return(list(M_stats = m_table, Occurrence_stats = o_table))
  })

  # returning final results
  names(var_stats) <- var_names
  return(var_stats)
}
