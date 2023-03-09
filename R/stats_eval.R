#' Statistics of environmental conditions in M and for occurrences (one variable)
#'
#' @description stats_eval helps in creating tables of descriptive statistics
#' of environmental conditions in accessible areas (M) and occurrence
#' records for one environmental variable at a time.
#'
#' @param stats (character) name or vector of names of functions to be applied
#' to get basic statistics of environmental values.
#' @param Ms a list of SpatVector objects representing the accessible area
#' (M) for each species to be analyzed. The order of species represented by each
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
#' @param variable a single SpatRaster layer of an environmental variable of
#' interest. See details.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 0.
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
#' @return
#' A list containing tables with statistics of the values in \code{variable},
#' for the species M and occurrences.
#'
#' @importFrom stats na.omit median
#' @importFrom terra extract crop
#'
#' @export
#'
#' @usage
#' stats_eval(stats = c("median", "range"), Ms, occurrences, species,
#'            longitude, latitude, variable, percentage_out = 0, verbose = TRUE)
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
#' # running stats
#' stat <- stats_eval(stats = c("mean", "sd", "median", "range", "quantile"),
#'                    Ms = m_list, occurrences = occ_list, species = "species",
#'                    longitude = "x", latitude = "y", variable = temp,
#'                    percentage_out = 0)

stats_eval <- function(stats = c("median", "range"), Ms, occurrences, species,
                       longitude, latitude, variable, percentage_out = 0,
                       verbose = TRUE) {
  # checking for potential errors
  if (missing(Ms)) {stop("Argument 'Ms' is missing.")}
  if (missing(occurrences)) {stop("Argument 'occurrences' is missing.")}
  if (missing(species)) {stop("Argument 'species' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(variable)) {stop("Argument 'variable' is missing.")}
  if (!is.list(Ms)) {stop("Argument 'Ms' must be a list.")}
  if (!is.list(occurrences)) {stop("Argument 'occurrences' must be a list.")}
  if (length(Ms) != length(occurrences)) {
    stop("'Ms' and 'occurrences' must have the same length and order of species listed must be the same.")
  }

  if (verbose == TRUE) {
    message("\nPreparing statistics from environmental layer and species data:")
  }

  sp_stats <- lapply(1:length(occurrences), function(j) {
    ## preparing e values
    mvar <- terra::crop(variable, Ms[[j]], mask = TRUE)
    mval <- na.omit(mvar[][, 1])
    if (percentage_out > 0) {
      medians <- median(mval)
      df_layer <- abs(mval - medians)
      names(df_layer) <- mval
      limit <- floor((100 - percentage_out) * length(df_layer)/100)
      df_layer <- sort(df_layer)[1:limit]
      mval <- as.numeric(names(df_layer))
    }

    occval <- na.omit(terra::extract(mvar,
                                     as.matrix(occurrences[[j]][, c(longitude,
                                                                    latitude)])))[, 1]

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

    if (verbose == TRUE) {
      message("\t", j, " of ", length(occurrences), " species finished")
    }
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
