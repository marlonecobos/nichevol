#' Bin table of environmental conditions in M and occurrences
#'
#' @description bin_table helps in creating a bin table of environmental
#' conditions in accessible areas (M) and species occurrence records
#' (i.e., table of characters).
#'
#' @param Ms a list of SpatialPolygons* objects representing the accessible area
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
#' @param variable a RasterLayer of an environmental variable of interest.
#' See details.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 5.
#' @param bin_size (numeric) size of bins. Range of environmental values to
#' be considered when creating each character in bin tables. See details.
#' Default = 10.
#'
#' @details
#' Coordinates in \code{occurrences}, SpatialPolygons* objects in \code{Ms}, and
#' RasterLayer in \code{variable} must coincide in the geographic projection in
#' which they are represented. WGS84 with no planar projection is recommended.
#'
#' Accessible area (M) is understood as the geographic area that has been
#' accessible for a species for relevant periods of time. Defining M is usually
#' a hard task, but also a very important one because it allows identifying
#' uncertainties about the ability of a species to maintain populations in
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) in \url{https://doi.org/10.1016/j.ecolmodel.2011.02.011}.
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
#' @importFrom raster extract crop mask
#'
#' @export
#'
#' @examples
#' # getting a variable at coarse resolution
#' \dontrun{
#' temp <- raster::getData("worldclim", var = "bio", res = 10)[[1]]
#'
#' # example data
#' data("m_list", package = "nichevol")
#' data("occ_list", package = "nichevol")
#'
#' # preparing bins
#' char_table <- bin_table(Ms = m_list, occurrences = occ_list, species = "species",
#'                         longitude = "x", latitude = "y", variable = temp,
#'                         percentage_out = 5, bin_size = 10)
#' }

bin_table <- function(Ms, occurrences, species, longitude, latitude, variable,
                      percentage_out = 5, bin_size = 10) {
  # checking for potential errors
  if (missing(Ms)) {stop("Argument Ms is missing.")}
  if (missing(occurrences)) {stop("Argument occurrences is missing.")}
  if (missing(species)) {stop("Argument species is missing.")}
  if (missing(longitude)) {stop("Argument longitude is missing.")}
  if (missing(latitude)) {stop("Argument latitude is missing.")}
  if (missing(variable)) {stop("Argument variable is missing.")}
  if (!is.list(Ms)) {stop("Argument Ms must be a list.")}
  if (!is.list(occurrences)) {stop("Argument occurrences must be a list.")}
  if (length(Ms) != length(occurrences)) {
    stop("Ms and occurrences must have the same length and order of species listed must be the same.")
  }

  M_range <- list()
  sp_range <- list()
  spnames <- vector()

  cat("\n   Preparing range values:\n")

  for (j in 1:length(occurrences)) {
    # processing
    ## get values of variable in M
    mvar <- raster::mask(raster::crop(variable, Ms[[j]]), Ms[[j]])
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
    occval <- na.omit(raster::extract(mvar, occurrences[[j]][, c(longitude,
                                                                 latitude)]))
    sp_range[[j]] <- range(occval)
    spnames[j] <- as.character(occurrences[[j]][1, species])

    cat("\t", j, "of", length(occurrences), "species finished\n")
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
  rownames(bin_table) <- gsub("_", " ", spnames)

  return(bin_table)
}
