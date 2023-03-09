#' Bin table of environmental conditions in M and for occurrences
#'
#' @description bin_table helps in creating a bin table of environmental
#' conditions in accessible areas (M) and for species occurrence records
#' (i.e., table of characters).
#'
#' @param Ms a list of SpatVector objects representing the accessible area
#' (M) for all species to be analyzed. The order of species represented by each
#' object here must coincide with the one in \code{occurrences}. See details.
#' @param occurrences a list of data.frames of occurrence records for all species.
#' The order of species represented by each data.frame must coincide with the one
#' in \code{Ms}. See details.
#' @param species (character) name of the column in occurrence data.frames that
#' contains the name of the species.
#' @param longitude (character) name of the column in occurrence files containing
#' values of longitude.
#' @param latitude (character) name of the column in occurrence files containing
#' values of latitude.
#' @param variable a single SpatRaster layer representing an environmental
#' variable of interest. See details.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 5.
#' @param n_bins (numeric) number of bins to be created from the range of
#' environmental values considered when creating each character in bin tables.
#' Default = 20. See details.
#' @param bin_size (numeric) argument deprecated, use n_bins instead.
#' @param verbose (logical) whether messages should be printed. Default = TRUE.
#'
#' @details
#' Coordinates in \code{occurrences}, SpatVector objects in \code{Ms}, and
#' SpatRaster in \code{variable} must coincide in the geographic projection in
#' which they are represented. WGS84 with no planar projection is recommended.
#'
#' Accessible area (M) is understood as the geographic area that has been
#' accessible for a species for relevant periods of time. Defining M is usually
#' a hard task, but also a very important one, because it allows identifying
#' uncertainties about the ability of a species to maintain populations in
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) <doi:10.1016/j.ecolmodel.2011.02.011>
#' and Machado-Stredel et al. (2021) <doi:10.21425/F5FBG48814>.
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
#' A list containing a table of characters to represent ecological niches of the
#' species of interest.
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
#' @importFrom stats na.omit median
#' @importFrom terra extract crop
#'
#' @export
#'
#' @usage
#' bin_table(Ms, occurrences, species, longitude, latitude, variable,
#'           percentage_out = 5, n_bins = 20, bin_size, verbose = TRUE)
#'
#' @examples
#' # example data
#' ## list of species records
#' data("occ_list", package = "nichevol")
#'
#' ## list of species accessible areas
#' m_files <- list.files(system.file("extdata", package = "nichevol"),
#'                       pattern = "m\\d.gpkg", full.names = TRUE)
#'
#' m_list <- lapply(m_files, terra::vect)
#'
#' ## raster variable
#' temp <- terra::rast(system.file("extdata", "temp.tif", package = "nichevol"))
#'
#'
#' # preparing bins
#' char_table <- bin_table(Ms = m_list, occurrences = occ_list, species = "species",
#'                         longitude = "x", latitude = "y", variable = temp,
#'                         percentage_out = 5, n_bins = 20)

bin_table <- function(Ms, occurrences, species, longitude, latitude, variable,
                      percentage_out = 5, n_bins = 20, bin_size, verbose = TRUE) {
  # checking for potential errors
  if (missing(Ms)) {stop("Argument 'Ms' is missing.")}
  if (missing(occurrences)) {stop("Argument 'occurrences' is missing.")}
  if (missing(species)) {stop("Argument 'species' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(variable)) {stop("Argument 'variable' is missing.")}
  if (!is.list(Ms)) {stop("Argument 'Ms' must be a list.")}
  if (!is.list(occurrences)) {stop("Argument 'occurrences' must be a list.")}
  if (!missing(bin_size)) {
    warning("Argument 'bin_size' is deprecated, using 'n_bins'.")
  }
  if (length(Ms) != length(occurrences)) {
    stop("'Ms' and 'occurrences' must have the same length and order of species listed must be the same.")
  }

  M_range <- list()
  sp_range <- list()
  spnames <- vector()

  if (verbose == TRUE) {
    message("\n   Preparing range values:")
  }

  for (j in 1:length(occurrences)) {
    # processing
    ## get values of variable in M
    mvar <- terra::crop(variable, Ms[[j]], mask = TRUE)
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
    occval <- na.omit(terra::extract(mvar,
                                     as.matrix(occurrences[[j]][, c(longitude,
                                                                    latitude)])))[, 1]
    sp_range[[j]] <- range(occval)
    spnames[j] <- as.character(occurrences[[j]][1, species])

    if (verbose == TRUE) {
      message("\t", j, " of ", length(occurrences), " species finished")
    }
  }

  # overall range
  M_range <- round(do.call(rbind, M_range), 2)
  sp_range <- round(do.call(rbind, sp_range), 2)
  overall_range <- range(c(c(M_range), c(sp_range)))

  # bin tables
  if (verbose == TRUE) {
    message("\n   Preparing bin tables using ranges:")
  }

  bin_size <- diff(overall_range) / n_bins

  bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
  rownames(bin_table) <- gsub("_", " ", spnames)

  return(bin_table)
}
