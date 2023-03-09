#' Statistics of environmental conditions in M and for occurrences (multiple variables)
#'
#' @description stats_evalues helps in creating csv files with statistics
#' of environmental conditions in accessible areas (M) and species occurrence
#' records. This is done using data read directly from a local directory, and
#' can be applied to various species and multiple variables.
#'
#' @param stats (character) name or vector of names of functions to be applied
#' to get basic statistics of environmental values.
#' @param M_folder (character) name of the folder containing files representing
#' the accessible area (M) for each species to be analyzed. See details.
#' @param M_format format of files representing the accessible area (M) for the
#' species. Names of M files must match the ones for occurrence files in
#' \code{occ_folder}. Format options are: "shp", "gpkg", or any of the options
#' supported by \code{\link[terra]{rast}} (e.g., "tif" or "asc").
#' @param occ_folder (character) name of the folder containing csv files of
#' occurrence data for all species. Names of csv files must match the ones of M
#' files in \code{M_folder}.
#' @param longitude (character) name of the column in occurrence files containing
#' values of longitude.
#' @param latitude (character) name of the column in occurrence files containing
#' values of latitude.
#' @param var_folder (character) name of the folder containing layers to
#' represent environmental variables.
#' @param var_format format of layers to represent environmental variables.
#' Format options are all the ones supported by \code{\link[terra]{rast}}
#' (e.g., "tif" or "asc").
#' @param round (logical) whether or not to round the values of one or more
#' variables after multiplying them times the value in \code{multiplication_factor}.
#' Default = FALSE. See details.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 0.
#' @param save (logical) whether or not to save the results in working directory.
#' Default = FALSE.
#' @param output_directory (character) name of the folder in which results will
#' be written.
#' @param overwrite (logical) whether or not to overwrite existing results in
#' \code{output_directory}. Default = FALSE.
#' @param verbose (logical) whether messages should be printed. Default = TRUE.
#'
#' @details
#' Coordinates in csv files in \code{occ_folder}, SpatVector-like files in
#' \code{M_folder}, and raster layers in \code{var_folder} must coincide in the
#' geographic projection in which they are represented. WGS84 with no planar
#' projection is recommended.
#'
#' Accessible area (M) is understood as the geographic area that has been
#' accessible for a species for relevant periods of time. Defining M is usually
#' a hard task, but also a very important one, because it allows identifying
#' uncertainties about the ability of a species to maintain populations in
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) <doi:10.1016/j.ecolmodel.2011.02.011>
#' and Machado-Stredel et al. (2021) <doi:10.21425/F5FBG48814>.
#'
#' Rounding variables may be useful when multiple variables are considered and
#' the values of some or all of them are too small (e.g., when using principal
#' components). To round specific variables arguments \code{round},
#' \code{round_names}, and \code{multiplication_factor}, must be used accordingly.
#'
#' The percentage to be defined in \code{percentage_out} excludes a percentage
#' of extreme environmental values to prevent the algorithm from considering
#' extremely rare environmental values in the accessible area for the species (M).
#' Being too rare, these values may have never been explored by the species;
#' therefore, including them in the process of preparation of the table of
#' characters (bin table) is risky.
#'
#' @return
#' A list named as the variables present in \code{var_folder}, containing all
#' tables with statistics of environmental values in M and in species records.
#' A folder named as in \code{output_directory} containing all resultant csv
#' files with the tables of statistics will be created if \code{save} is set as
#' TRUE.
#'
#' @importFrom stats na.omit median
#' @importFrom utils write.csv read.csv
#' @importFrom terra extract crop nlyr rast vect
#'
#' @export
#'
#' @usage
#' stats_evalues(stats = c("median", "range"), M_folder, M_format, occ_folder,
#'               longitude, latitude, var_folder, var_format, round = FALSE,
#'               round_names, multiplication_factor = 1, percentage_out = 0,
#'               save = FALSE, output_directory, overwrite = FALSE,
#'               verbose = TRUE)
#'
#' @examples
#' # preparing data and directories for examples
#' ## directories
#' tempdir <- file.path(tempdir(), "nevol_test")
#' dir.create(tempdir)
#'
#' cvariables <- paste0(tempdir, "/variables")
#' dir.create(cvariables)
#'
#' records <- paste0(tempdir, "/records")
#' dir.create(records)
#'
#' m_areas <- paste0(tempdir, "/M_areas")
#' dir.create(m_areas)
#'
#' ## data
#' data("occ_list", package = "nichevol")
#'
#' temp <- system.file("extdata", "temp.tif", package = "nichevol")
#'
#' m_files <- list.files(system.file("extdata", package = "nichevol"),
#'                       pattern = "m\\d.gpkg", full.names = TRUE)
#'
#' ## writing data in temporal directories
#' spnames <- sapply(occ_list, function (x) as.character(x[1, 1]))
#' ocnames <-  paste0(records, "/", spnames, ".csv")
#'
#' occs <- lapply(1:length(spnames), function (x) {
#'   write.csv(occ_list[[x]], ocnames[x], row.names = FALSE)
#' })
#'
#' to_replace <- paste0(system.file("extdata", package = "nichevol"), "/")
#'
#' otemp <- gsub(to_replace, "", temp)
#' file.copy(from = temp, to = paste0(cvariables, "/", otemp))
#'
#' file.copy(from = m_files, to = paste0(m_areas, "/", spnames, ".gpkg"))
#' stats <- stats_evalues(stats = c("median", "range"), M_folder = m_areas,
#'                        M_format = "gpkg", occ_folder = records,
#'                        longitude = "x", latitude = "y",
#'                        var_folder = cvariables, var_format = "tif",
#'                        percentage_out = 5)

stats_evalues <- function(stats = c("median", "range"), M_folder, M_format,
                          occ_folder, longitude, latitude, var_folder,
                          var_format, round = FALSE, round_names,
                          multiplication_factor = 1, percentage_out = 0,
                          save = FALSE, output_directory, overwrite = FALSE,
                          verbose = TRUE) {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
  if (save == TRUE) {
    if (missing(output_directory)) {
      stop("Argument 'output_directory' is missing.")
    } else {
      if (overwrite == FALSE & dir.exists(output_directory)) {
        stop("'output_directory' already exists, to replace it use 'overwrite' = TRUE.")
      }
      if (overwrite == TRUE & dir.exists(output_directory)) {
        unlink(x = output_directory, recursive = TRUE, force = TRUE)
      }
    }
  }

  # formats and data to start
  if (verbose == TRUE) {
    message("\nPreparing data, please wait...\n")
  }

  ## records
  occlist <- list.files(path = occ_folder, pattern = ".csv$", full.names = TRUE)

  ## M areas
  M_patt <- paste0(M_format, "$")
  mlist <- list.files(path = M_folder, pattern = M_patt, full.names = TRUE)

  ## species names
  subs <- paste0(".", M_format)
  spnames <- gsub(subs, "", list.files(path = M_folder, pattern = M_patt))

  ## variables
  v_patt <- paste0(var_format, "$")
  variables <- terra::rast(list.files(path = var_folder, pattern = v_patt,
                                      full.names = TRUE))

  var_names <- names(variables)
  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    noround <- var_names[!var_names %in% round_names]
    variables <- c(variables[[noround]], rounds)
  }

  # directory for results
  if (save == TRUE) {dir.create(output_directory)}
  if (verbose == TRUE) {
    message("Preparing statistics from environmental layers and species data:")
  }

  n_vars <- terra::nlyr(variables)

  var_stats <- lapply(1:n_vars, function(i) {
    sp_stats <- lapply(1:length(occlist), function(j) {
      ## M
      if (M_format %in% c("shp", "gpkg")) {
        M <- terra::vect(mlist[j])
      } else {
        M <- terra::rast(mlist[j])
      }

      ## occurrences
      occ <- read.csv(occlist[j])

      # processing
      ## get values of variables in M
      mvar <- terra::crop(variables[[i]], M, mask = TRUE)
      mval <- na.omit(mvar[][, 1])

      if (percentage_out > 0) {
        medians <- median(mval)
        df_layer <- abs(mval - medians)
        names(df_layer) <- mval
        limit <- floor((100 - percentage_out) * length(df_layer) / 100)
        df_layer <- sort(df_layer)[1:limit]
        mval <- as.numeric(names(df_layer))
      }

      ## occurrences
      occval <- na.omit(terra::extract(mvar, as.matrix(occ[, c(longitude,
                                                               latitude)]))[, 1])

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

      message("\t", j, " of ", length(occlist), " species finished")
      return(list(M = unlist(m_stats), Occurrences = unlist(o_stats)))
    })

    # preparing tables with results
    spnames <- gsub("_", " ", spnames)
    m_table <- data.frame(Species = spnames,
                          do.call(rbind, lapply(sp_stats, function(x) {x[[1]]})))
    o_table <- data.frame(Species = spnames,
                          do.call(rbind, lapply(sp_stats, function(x) {x[[2]]})))

    # write table
    if (save == TRUE) {
      write.csv(m_table,
                paste0(output_directory, "/", var_names[i], "_M_stats.csv"),
                row.names = FALSE)
      write.csv(o_table,
                paste0(output_directory, "/", var_names[i], "_Occurrence_stats.csv"),
                row.names = FALSE)
    }

    if (verbose == TRUE) {
      message(i, " of ", n_vars, " variables processed")
    }

    return(list(M_stats = m_table, Occurrence_stats = o_table))
  })

  # returning final results
  names(var_stats) <- var_names
  return(var_stats)
}
