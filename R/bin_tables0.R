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
#' @param round (logical) whether or not to round the values of one or more
#' variables after multiplying them times the value in \code{multiplication_factor}.
#' Default = FALSE. See details.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 5.
#' @param bin_size (numeric) size of bins. Range of environmental values to
#' be considered when creating each character in bin tables. See details.
#' Default = 10.
#' @param save (logical) whether or not to save the results in working directory.
#' Default = FALSE.
#' @param overwrite (logical) whether or not to overwrite existing results in
#' \code{output_directory}. Default = FALSE.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Species_E_bins".
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
#' uncertainties about the ability of a species to maintain populations in
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) in \url{https://doi.org/10.1016/j.ecolmodel.2011.02.011}.
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
#' The argument \code{bin_size} helps to create characters that represent not
#' only one value of an environmental variable, but a range of environmental
#' conditions. For instance, if a variable of precipitation in mm is used, a
#' value of 10 for \code{bin_size} indicates that each character will represent
#' a class that correspond to 10 continuous values of precipitation (e.g., from
#' 100 to 110 mm).
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
#' @importFrom raster extract crop mask nlayers raster stack
#' @importFrom rgdal readOGR
#'
#' @export
#'
#' @usage
#' bin_tables0(M_folder, M_format, occ_folder, longitude,
#'   latitude, var_folder, var_format, round = FALSE,
#'   round_names, multiplication_factor = 1,
#'   percentage_out = 5, bin_size = 10, save = FALSE,
#'   overwrite = FALSE, output_directory = "Species_E_bins")
#'
#' @examples
#' # example of how to define arguments, check argument descriptions above
#' \dontrun{
#' bins <- bin_tables0(M_folder = "Folder_with_Ms", M_format = "shp",
#'                     occ_folder = "Folder_with_occs", longitude = "lon_column",
#'                     latitude = "lat_column", var_folder = "Folder_with_vars",
#'                     var_format = "GTiff", percentage_out = 5, bin_size = 10)
#' }

bin_tables0 <- function(M_folder, M_format, occ_folder, longitude,
                        latitude, var_folder, var_format,
                        round = FALSE, round_names, multiplication_factor = 1,
                        percentage_out = 5, bin_size = 10, save = FALSE,
                        overwrite = FALSE, output_directory = "Species_E_bins") {
  # checking for potential errors
  if (missing(M_folder)) {stop("Argument 'M_folder' is missing.")}
  if (missing(M_format)) {stop("Argument 'M_format' is missing.")}
  if (missing(occ_folder)) {stop("Argument 'occ_folder' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
  if (save == TRUE) {
    if (overwrite == FALSE & dir.exists(output_directory)) {
      stop("'output_directory' already exists, to replace it use overwrite = TRUE.")
    }
    if (overwrite == TRUE & dir.exists(output_directory)) {
      unlink(x = output_directory, recursive = TRUE, force = TRUE)
    }
  }

  # formats and data to start
  message("\nPreparing data, please wait...\n\n")
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
  if (save == TRUE) {dir.create(output_directory)}
  message("Preparing range values and bin tables from environmental layers and species data:\n")
  nvars <- raster::nlayers(variables)

  bin_tabs <- lapply(1:nvars, function(i) {
    # data
    M_range <- list()
    sp_range <- list()

    message("\n   Preparing range values:\n")

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

      message("\t", j, " of ", length(occlist), " species finished\n")
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
    message("   Preparing bin tables using ranges:\n")

    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    rownames(bin_table) <- gsub("_", " ", spnames)

    # write table
    if (save == TRUE) {
      write.csv(bin_table,
                paste0(output_directory, "/", names(variables)[i], "_bin_table.csv"),
                row.names = TRUE)
    }

    message(i, " of ", nvars, " variables processed\n")
    return(bin_table)
  })

  names(bin_tabs) <- names(variables)
  return(bin_tabs)
}
