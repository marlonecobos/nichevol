#' Histograms of environmental conditions in M and for occurrences
#'
#' @description histograms_env creates PDF files with histogram plots of
#' environmental conditions in M, lines for the confidence limits of values in
#' M, and the location of values in occurrence records. This is done using data
#' read directly from a local directory, and can be applied to various species
#' and multiple variables.
#'
#' @param M_folder (character) name of the folder containing files representing
#' the accessible area (M) for all species to be analyzed. See details.
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
#' @param CL_lines (numeric) confidence limits of environmental values in M to
#' be plotted as lines in the histograms. See details. Default = c(95, 99).
#' @param col colors for lines representing confidence limits. If NULL, colors
#' are selected from a gray palette. Default = NULL.
#' @param round (logical) whether or not to round values of one or more
#' variables after multiplying them times the value in \code{multiplication_factor}.
#' Default = FALSE. See details.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param save_ranges (logical) whether or not to save the values identified as
#' ranges considering the whole set of values and confidence limits defined in
#' \code{CL_lines}. Default = FALSE.
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
#' uncertainties about the ability of a species to maintain populations under
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) <doi:10.1016/j.ecolmodel.2011.02.011>
#' and Machado-Stredel et al. (2021) <doi:10.21425/F5FBG48814>.
#'
#' Rounding variables may be useful when multiple variables are considered and
#' the values of some or all of them are too small (e.g., when using principal
#' components). To round specific variables arguments \code{round},
#' \code{round_names}, and \code{multiplication_factor}, must be used accordingly.
#'
#' @return
#' A list of data.frames containing intervals of environmental values in species
#' occurrences and accessible areas (M), as well as values corresponding to the
#' confidence limits defined in \code{CL_lines}. A folder named as
#' in \code{output_directory} containing all resulting PDF files (one per
#' variable) with histograms for all species. Files (csv) of ranges found during
#' the analyses will be also written in \code{output_directory} if
#' \code{save_ranges} is set as TRUE.
#'
#' @importFrom grDevices gray.colors
#' @importFrom utils write.csv read.csv
#' @importFrom stats na.omit median
#' @importFrom terra extract crop nlyr rast vect
#'
#' @export
#'
#' @usage
#' histograms_env(M_folder, M_format, occ_folder, longitude, latitude,
#'                var_folder, var_format, CL_lines = c(95, 99), col = NULL,
#'                round = FALSE, round_names = NULL, multiplication_factor = 1,
#'                save_ranges = FALSE, output_directory, overwrite = FALSE,
#'                verbose = TRUE)
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
#' histdir <- paste0(tempdir, "/Hists")
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
#'
#' # running analysis to produce plots
#' hists <- histograms_env(M_folder = m_areas, M_format = "gpkg",
#'                         occ_folder = records, longitude = "x",
#'                         latitude = "y", var_folder = cvariables,
#'                         var_format = "tif", output_directory = histdir)

histograms_env <- function(M_folder, M_format, occ_folder, longitude,
                           latitude, var_folder, var_format,
                           CL_lines = c(95, 99), col = NULL, round = FALSE,
                           round_names = NULL, multiplication_factor = 1,
                           save_ranges = FALSE, output_directory,
                           overwrite = FALSE, verbose = TRUE) {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
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

  lcll <- length(CL_lines)
  if (is.null(col)) {
    col <- sort(gray.colors(lcll + 1), decreasing = TRUE)[1:lcll]
  }
  if (round == TRUE & is.null(round_names)) {
    stop("Argument 'round_names' cannot be NULL if round = TRUE.")
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

  ### rounding variables
  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- c(variables[[noround]], rounds)
  }

  # directory for results
  dir.create(output_directory)

  nvars <- terra::nlyr(variables)

  if (verbose == TRUE) {
    message("Preparing environmental values and histograms from layers and species data:")
  }

  ranges <- lapply(1:nvars, function(i) {
    # data
    df_layer <- list()
    occ_dfs <- list()

    M_ranges <- list()
    M_limits <- list()
    sp_ranges <- list()
    y_values <- list()

    if (verbose == TRUE) {
      message("\n   Preparing environmental values:")
    }

    for (j in 1:length(occlist)) {
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

      ## distance of each absolute value to median value
      medians <- median(mval)
      df_layer[[j]] <- abs(mval - medians)
      names(df_layer[[j]]) <- mval

      occval <- na.omit(terra::extract(mvar,
                                       as.matrix(occ[, c(longitude, latitude)])))[, 1]
      occ_dfs[[j]] <- abs(occval - medians)
      names(occ_dfs[[j]]) <- occval

      y_values[[j]] <- c(max(table(mval)), max(table(df_layer[[j]])))

      ## ranges of real values
      M_limit <- lapply(1:lcll, function(k) {
        limit <- floor((CL_lines[k]) * length(df_layer[[j]]) / 100)
        df_layera <- sort(df_layer[[j]])[1:limit]

        range(as.numeric(names(df_layera)))
      })

      M_ranges[[j]] <- range(mval)
      M_limits[[j]] <- do.call(c, M_limit)
      sp_ranges[[j]] <- range(occval)

      if (verbose == TRUE) {
        message("\t", j, " of ", length(occlist), " species finished")
      }
    }

    ## ranges final values for variables
    limits <- do.call(rbind, M_limits)
    ranges <- data.frame(gsub("_", " ", spnames), do.call(rbind, sp_ranges),
                         do.call(rbind, M_ranges), limits)
    colnames(ranges) <- c("Species", "Species_lower", "Species_upper",
                               "M_lower", "M_upper", paste0("M_",
                                                            rep(CL_lines, each = 2),
                                                            c("_lowerCL", "_upperCL")))

    if (save_ranges == TRUE) {
      write.csv(ranges, file = paste0(output_directory, "/Ranges_",
                                      names(variables)[i], ".csv"),
                row.names = FALSE)
    }

    ## frecuency of each value in M
    if (verbose == TRUE) {
      message("\n   Preparing histogram plots using environmental values...")
    }
    pdf_histograms(env_data = df_layer, occ_data = occ_dfs, y_values = y_values,
                   sp_names = spnames, variable_name = names(variables)[i],
                   CL_lines = CL_lines, limits = limits, col = col,
                   output_directory = output_directory)

    message(i, " of ", nvars, " variables processed")
    return(ranges)
  })

  if (save_ranges == TRUE) {
    if (verbose == TRUE) {
      message("\ncsv files with the environmental ranges were saved in ",
              output_directory, "\n")
    }
  }

  names(ranges) <- names(variables)

  return(ranges)
}
