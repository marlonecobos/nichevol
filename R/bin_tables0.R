#' Bin tables of environmental conditions in M and for occurrences from data
#'
#' @description bin_tables0 helps in creating bin tables of environmental
#' conditions in accessible areas (M) and species occurrence records
#' (i.e., table of characters). This is done using data read directly from a
#' local directory, and can be applied to various species and multiple variables.
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
#' @param round (logical) whether or not to round the values of one or more
#' variables after multiplying them times the value in \code{multiplication_factor}.
#' Default = FALSE. See details.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 5.
#' @param n_bins (numeric) number of bins to be created from the range of
#' environmental values considered when creating each character in bin tables.
#' Default = 20. See details.
#' @param bin_size (numeric) argument deprecated, use n_bins instead.
#' @param save (logical) whether or not to save the results in working directory.
#' Default = FALSE.
#' @param output_directory (character) name of the folder in which results will be
#' written.
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
#' of extreme environmental values to prevent from considering extremely rare
#' environmental values in the accessible area for the species (M). Being too
#' rare, these values may have never been explored by the species; therefore,
#' including them in the process of preparation of the table of characters
#' (bin table) is risky.
#'
#' The argument \code{n_bins} helps to define how many characters (bins) will be
#' considered for the range of values in each variable. This is, a value of 20
#' determines that a range of temperature (5-25) will be split approximately
#' every 1 degree. The argument \code{bin_size} has been deprecated.
#'
#' @return
#' A list named as the variables present in \code{var_folder}, containing all
#' tables of characters. A folder named as in \code{output_directory} containing
#' all resultant csv files with the tables of characters will be created if
#' \code{save} is set as TRUE.
#'
#' Potential values for characters are:
#' - "1" = the species is present in those environmental conditions.
#' - "0" = the species is not present in those environmental conditions. This is,
#' those environmental conditions inside the accessible area (M) are more extreme
#' than the ones used for the species.
#' - "?" = there is no certainty about the species presence in those environmental
#' conditions. This happens in environmental combinations more extreme than the
#' ones found in the accessible area (M), when environmental conditions in
#' species records are as extreme as the most extreme ones in M.
#'
#' @importFrom utils write.csv read.csv
#' @importFrom stats na.omit median
#' @importFrom terra extract crop nlyr rast vect
#'
#' @export
#'
#' @usage
#' bin_tables0(M_folder, M_format, occ_folder, longitude,
#'             latitude, var_folder, var_format, round = FALSE,
#'             round_names, multiplication_factor = 1,
#'             percentage_out = 5, n_bins = 20, bin_size, save = FALSE,
#'             output_directory, overwrite = FALSE, verbose = TRUE)
#'
#' @examples
#' # preparing data and directories for example
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
#'
#' # preparing tables
#' tabs <- bin_tables0(M_folder = m_areas, M_format = "gpkg", occ_folder = records,
#'                     longitude = "x", latitude = "y", var_folder = cvariables,
#'                     var_format = "tif")


bin_tables0 <- function(M_folder, M_format, occ_folder, longitude,
                        latitude, var_folder, var_format,
                        round = FALSE, round_names, multiplication_factor = 1,
                        percentage_out = 5, n_bins = 20, bin_size, save = FALSE,
                        output_directory, overwrite = FALSE, verbose = TRUE) {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
  if (!missing(bin_size)) {
    warning("Argument 'bin_size' is deprecated, using 'n_bins'.")
  }
  if (save == TRUE) {
    if (missing(output_directory)) {
      stop("Argument 'output_directory' is missing.")
    } else {
      if (overwrite == FALSE & dir.exists(output_directory)) {
        stop("'output_directory' already exists, to replace it use overwrite = TRUE.")
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

  ### rounding variables
  if (round == TRUE) {
    rounds <- round(variables[[round_names]] * multiplication_factor)
    var_names <- names(variables)
    noround <- var_names[!var_names %in% round_names]
    variables <- c(variables[[noround]], rounds)
  }

  # directory for results
  if (save == TRUE) {dir.create(output_directory)}
  if (verbose == TRUE) {
    message("Preparing range values and bin tables from environmental layers and species data:")
  }

  nvars <- terra::nlyr(variables)

  bin_tabs <- lapply(1:nvars, function(i) {
    # data
    M_range <- list()
    sp_range <- list()

    if (verbose == TRUE) {
      message("\n   Preparing range values:")
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
      df_layer <- abs(mval - medians)
      names(df_layer) <- mval

      ## limit
      limit <- floor((100 - percentage_out) * length(df_layer) / 100)
      df_layer <- sort(df_layer)[1:limit]

      M_range[[j]] <- range(as.numeric(names(df_layer)))

      ## occurrences
      occval <- na.omit(terra::extract(mvar, as.matrix(occ[, c(longitude,
                                                      latitude)])))[, 1]
      sp_range[[j]] <- range(occval)

      if (verbose == TRUE) {
        message("\t", j, " of ", length(occlist), " species finished")
      }
    }

    # overall range
    M_range <- round(do.call(rbind, M_range), 2)
    sp_range <- round(do.call(rbind, sp_range), 2)
    overall_range <- range(c(c(M_range), c(sp_range)))

    # bin tables
    if (verbose == TRUE) {
      message("\n   Preparing character tables using ranges:")
    }

    bin_size <- diff(overall_range) / n_bins

    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    rownames(bin_table) <- gsub("_", " ", spnames)

    # write table
    if (save == TRUE) {
      write.csv(bin_table,
                paste0(output_directory, "/", names(variables)[i], "_bin_table.csv"),
                row.names = TRUE)
    }

    if (verbose == TRUE) {
      message(i, " of ", nvars, " variables processed")
    }
    return(bin_table)
  })

  names(bin_tabs) <- names(variables)
  return(bin_tabs)
}
