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
#' in \code{\link[raster]{writeFormats}} (e.g., "GTiff").
#' @param occ_folder (character) name of the folder containing csv files of
#' occurrence data for all species. Names of csv files must match the ones of M
#' files in \code{M_folder}.
#' @param longitude (character) name of the column in occurrence files containing
#' values of longitude.
#' @param latitude (character) name of the column in occurrence files containing
#' values of latitude.
#' @param var_folder (character) name of the folder containing layers to
#' represent environmental variables.
#' @param var_format format of layers to represent environmental variables. See
#' options in \code{\link[raster]{writeFormats}} (e.g., "GTiff").
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
#' @param overwrite (logical) whether or not to overwrite existing results in
#' \code{output_directory}. Default = FALSE.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Histogram_ranges".
#'
#' @details
#' Coordinates in csv files in \code{occ_folder}, SpatialPolygons*-like files in
#' \code{M_folder}, and raster layers in \code{var_folder} must coincide in the
#' geographic projection in which they are represented. WGS84 with no planar
#' projection is recommended.
#'
#' Accessible area (M) is understood as the geographic area that has been
#' accessible for a species for relevant periods of time. Defining M is usually
#' a hard task, but also a very important one, because it allows identifying
#' uncertainties about the ability of a species to maintain populations under
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) in \url{https://doi.org/10.1016/j.ecolmodel.2011.02.011}.
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
#' @importFrom raster extract crop mask nlayers raster stack
#' @importFrom rgdal readOGR
#'
#' @export
#'
#' @usage
#' histograms_env(M_folder, M_format, occ_folder, longitude, latitude,
#'   var_folder, var_format, CL_lines = c(95, 99), col = NULL,
#'   round = FALSE, round_names = NULL, multiplication_factor = 1,
#'   save_ranges = FALSE, overwrite = FALSE,
#'   output_directory = "Histogram_ranges")
#'
#' @examples
#' # example of how to define arguments, check argument descriptions above
#' \dontrun{
#' hists <- histograms_env(M_folder = "Folder_with_Ms", M_format = "shp",
#'                         occ_folder = "Folder_with_occs", longitude = "lon_column",
#'                         latitude = "lat_column", var_folder = "Folder_with_vars",
#'                         var_format = "GTiff")
#' }

histograms_env <- function(M_folder, M_format, occ_folder, longitude,
                           latitude, var_folder, var_format,
                           CL_lines = c(95, 99), col = NULL, round = FALSE,
                           round_names = NULL, multiplication_factor = 1,
                           save_ranges = FALSE, overwrite = FALSE,
                           output_directory = "Histogram_ranges") {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
  if (overwrite == FALSE & dir.exists(output_directory)) {
    stop("'output_directory' already exists, to replace it use overwrite = TRUE.")
  }
  if (overwrite == TRUE & dir.exists(output_directory)) {
    unlink(x = output_directory, recursive = TRUE, force = TRUE)
  }

  lcll <- length(CL_lines)
  if (is.null(col)) {
    col <- sort(gray.colors(lcll + 1), decreasing = TRUE)[1:lcll]
  }
  if (round == TRUE & is.null(round_names)) {
    stop("Argument 'round_names' cannot be NULL if round = TRUE.")
  }

  # formats and data to start
  message("\nPreparing data, please wait...\n")
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

  nvars <- raster::nlayers(variables)
  message("Preparing environmental values and histograms from layers and species data:")

  ranges <- lapply(1:nvars, function(i) {
    # data
    df_layer <- list()
    occ_dfs <- list()

    M_ranges <- list()
    M_limits <- list()
    sp_ranges <- list()
    y_values <- list()

    message("\n   Preparing environmental values:")

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
      df_layer[[j]] <- abs(mval - medians)
      names(df_layer[[j]]) <- mval

      occval <- na.omit(raster::extract(mvar, occ[, c(longitude, latitude)]))
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

      message("\t", j, " of ", length(occlist), " species finished")
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
    message("\n   Preparing histogram plots using environmental values...")
    pdf_histograms(env_data = df_layer, occ_data = occ_dfs, y_values = y_values,
                   sp_names = spnames, variable_name = names(variables)[i],
                   CL_lines = CL_lines, limits = limits, col = col,
                   output_directory = output_directory)

    message(i, " of ", nvars, " variables processed")
    return(ranges)
  })

  if (save_ranges == TRUE) {
    message("\ncsv files with the environmental ranges were saved in ",
        output_directory, "\n")
  }

  names(ranges) <- names(variables)

  return(ranges)
}
