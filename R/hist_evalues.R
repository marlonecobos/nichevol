#' Histograms of environmental conditions in M and for occurrences (one species)
#'
#' @description hist_evalues helps in creating histograms to explore environmental
#' conditions in M, lines for the confidence limits of values in M, and the
#' location of values in occurrence records, for one species at the time.
#'
#' @param M a SpatVector object representing the accessible area (M)
#' for one species. See details.
#' @param occurrences a data.frame of occurrence records for one species. See
#' details.
#' @param species (character) name of the column in \code{occurrences} that
#' contains the name of the species.
#' @param longitude (character) name of the column in \code{occurrences} containing
#' values of longitude.
#' @param latitude (character) name of the column in \code{occurrences} containing
#' values of latitude.
#' @param variable a single SpatRaster layer representing an environmental
#' variable of interest. See details.
#' @param CL_lines (numeric) confidence limits of environmental values in M to
#' be plotted as lines in the histograms. See details. Default = c(95, 99).
#' @param col colors for lines representing confidence limits. If NULL, colors
#' are selected from a gray palette. Default = NULL.
#'
#' @details
#' Coordinates in \code{occurrences}, SpatVector object in \code{M}, and
#' SpatRaster in \code{variable} must coincide in the geographic projection in
#' which they are represented. WGS84 with no planar projection is recommended.
#'
#' The accessible area (M) is understood as the geographic area that has been
#' accessible to a species over relevant periods of time. Defining M is usually
#' a hard task, but also a very important one because it allows identifying
#' uncertainties about the ability of a species to maintain populations under
#' certain environmental conditions. For further details on this topic, see
#' Barve et al. (2011) <doi:10.1016/j.ecolmodel.2011.02.011>
#' and Machado-Stredel et al. (2021) <doi:10.21425/F5FBG48814>.
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics abline hist points
#' @importFrom stats na.omit median
#' @importFrom terra extract crop
#'
#' @export
#'
#' @usage
#' hist_evalues(M, occurrences, species, longitude, latitude, variable,
#'              CL_lines = c(95, 99), col = NULL)
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
#' hist_evalues(M = m_list[[1]], occurrences = occ_list[[1]], species = "species",
#'              longitude = "x", latitude = "y", variable = temp,
#'              CL_lines = c(95, 99), col = c("blue", "red"))

hist_evalues <- function(M, occurrences, species, longitude, latitude, variable,
                         CL_lines = c(95, 99), col = NULL) {
  # checking for potential errors
  if (missing(M)) {stop("Argument 'M' is missing.")}
  if (missing(occurrences)) {stop("Argument 'occurrences' is missing.")}
  if (missing(species)) {stop("Argument 'species' is missing.")}
  if (missing(longitude)) {stop("Argument 'longitude' is missing.")}
  if (missing(latitude)) {stop("Argument 'latitude' is missing.")}
  if (missing(variable)) {stop("Argument 'variable' is missing.")}
  if (is.null(col)) {
    col <- sort(gray.colors(lcll + 1), decreasing = TRUE)[1:lcll]
  }
  col1 <- rep(col, each = 2)

  # par settings
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # preparing data for plotting
  ## species name
  spname <- as.character(occurrences[1, species])

  ## confidence limits length
  lcll <- length(CL_lines)

  ## get values of variables in M
  mvar <- terra::crop(variable, M, mask = TRUE)
  mval <- na.omit(mvar[][, 1])

  ## distance of each absolute value to median value
  medians <- median(mval)
  df_layer <- abs(mval - medians)
  names(df_layer) <- mval

  occval <- na.omit(terra::extract(mvar,
                                   as.matrix(occurrences[, c(longitude,
                                                             latitude)])))[, 1]
  occ_dfs <- abs(occval - medians)
  names(occ_dfs) <- occval

  y_values <- c(max(table(mval)), max(table(df_layer)))

  ## ranges of real values
  M_limit <- lapply(1:lcll, function(k) {
    limit <- floor((CL_lines[k]) * length(df_layer) / 100)
    df_layera <- sort(df_layer)[1:limit]

    range(as.numeric(names(df_layera)))
  })

  M_ranges <- range(mval)
  M_limits <- do.call(c, M_limit)
  sp_ranges <- range(occval)

  ## for legend
  sym_legend <- c("", "", "Occurrences", paste0(CL_lines, "% CI"))
  lin <- c(NA, NA, NA, rep(1, length(CL_lines)))
  poi <- c(NA, NA, 1, rep(NA, length(CL_lines)))
  colss <- c(NA, NA, "darkgreen", col)

  # plotting
  par(mfrow = c(1, 2), cex = 0.8, mar = (c(4.5, 4, 3.5, 1) + 0.1))
  ## actual values
  hist(as.numeric(names(df_layer)),
       main = gsub("_", " ", spname), xlab = "Variable values")
  points(as.numeric(names(occ_dfs)),
         rep(y_values[1], length(occ_dfs)), col = "darkgreen",
         pch = 1)
  abline(v = M_limits, col = col1)

  ## median deviation
  hist(df_layer, main = gsub("_", " ", spname),
       xlab = "Median deviation")
  points(occ_dfs, rep(y_values[2], length(occ_dfs)),
         col = "darkgreen", pch = 1)
  limit <- floor(CL_lines * length(df_layer) / 100)
  abline(v = sort(df_layer)[limit], col = col)

  ## legend
  legend("topright", legend = sym_legend, lty = lin, pch = poi, col = colss,
         cex = 0.9, box.col = "white", bg = "white", inset = 0)
}
