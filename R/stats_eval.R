#' Statistics of environmental conditions in M and occurrences (one variable)
#'
#' @description stats_eval helps in creating tables of descriptive statistics
#' of environmental conditions in accessible areas (M) and species occurrence
#' records for one environmental variable at the time.
#'
#' @param stats (character) name or vector of names of functions to be applied
#' to get basic statistics of environmental values.
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
#' to be excluded in bin creation for further analyses. See details. Default = 0.
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
#' @return
#' A list containing tables with statistics of the values in \code{variable},
#' for the species M and occurrences.
#'
#' @importFrom stats na.omit median
#' @importFrom raster extract crop mask
#'
#' @export
#'
#' @examples
#' # getting a variable at coarse resolution
#' temp <- getData("worldclim", var = "bio", res = 10)[[1]]
#'
#' # example data
#' data("m_list", package = "nichevol")
#' data("occ_list", package = "nichevol")
#'
#' # running stats
#' stat <- stats_eval(stats = c("mean", "sd", "median", "range", "quantile"),
#'                    Ms = m_list, occurrences = occ_list, species = "species",
#'                    longitude = "x", latitude = "y", variable = temp,
#'                    percentage_out = 0)

stats_eval <- function(stats = c("median", "range"), Ms, occurrences, species,
                       longitude, latitude, variable, percentage_out = 0) {
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

  cat("Preparing statistics from environmental layer and species data:\n")
  sp_stats <- lapply(1:length(occurrences), function(j) {
    ## preparing e values
    mvar <- raster::mask(raster::crop(variable, Ms[[j]]), Ms[[j]])
    mval <- na.omit(mvar[])
    if (percentage_out > 0) {
      medians <- median(mval)
      df_layer <- abs(mval - medians)
      names(df_layer) <- mval
      limit <- floor((100 - percentage_out) * length(df_layer)/100)
      df_layer <- sort(df_layer)[1:limit]
      mval <- as.numeric(names(df_layer))
    }

    occval <- na.omit(raster::extract(mvar, occurrences[[j]][, c(longitude,
                                                                 latitude)]))

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

    spn <- as.character(occurrences[[j]][1, species])

    cat("\t", j, "of", length(occurrences), "species finished\n")
    return(list(sp = spn, M = unlist(m_stats), Occurrences = unlist(o_stats)))
  })

  # preparing tables with results
  spnames <- gsub("_", " ", unlist(lapply(sp_stats, function(x) {x[[1]]})))
  m_table <- data.frame(Species = spnames,
                        do.call(rbind, lapply(sp_stats, function(x) {x[[2]]})))
  o_table <- data.frame(Species = spnames,
                        do.call(rbind, lapply(sp_stats, function(x) {x[[3]]})))

  return(list(M_stats = m_table, Occurrence_stats = o_table))
}
