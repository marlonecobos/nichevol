#' Bin tables of environmental conditions for virtual species
#'
#' @description bin_tables_virtual helps create csv files of environmental
#' condition bin tables for virtual species. It can do this for multiple
#' species and multiple variables.
#'
#' @param virtualsp_folder (character) name of the folder containing files
#' representing the virtually occupied area of all species to be analyzed.
#' @param virtualsp_format format of files representing the virtually occupied
#' area of all species. Options are: "shp", "gpkg", or any of the options in
#' \code{\link[raster]{writeFormats}}.
#' @param var_folder (character) name of the folder containing layers to
#' represent environmental variables.
#' @param var_format format of layers to represent environmental variables. See
#' options in \code{\link[raster]{writeFormats}}.
#' @param round (logical) whether or not to round the values of one or more
#' variables after multiplying them times \code{multiplication_factor}.
#' Default = FALSE.
#' @param round_names (character) names of the variables to be rounded.
#' Default = NULL. If \code{round} = TRUE, names must be defined.
#' @param multiplication_factor (numeric) value to be used to multiply the
#' variables defined in \code{round_names}. Default = 1.
#' @param bin_size (numeric) size of bins. Interval of values to be considered
#' when creating bin tables. Default = 10.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Virtual_species_E_bins".
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all resultant csv
#' files for all variables, with bins for all species. Results will also be
#' returned as a list.

bin_tables_virtual <- function(virtualsp_folder, virtualsp_format, var_folder, var_format,
                            round = FALSE, round_names = NULL,
                            multiplication_factor = 1, bin_size = 10,
                            output_directory = "Virtual_species_E_bins"){
  # checking for potential errors
  if (missing(virtualsp_folder)) {stop("Argument 'virtualsp_folder' is missing.")}
  if (missing(virtualsp_format)) {stop("Argument 'virtualsp_format' is missing.")}
  if (missing(var_folder)) {stop("Argument 'var_folder' is missing.")}
  if (missing(var_format)) {stop("Argument 'var_format' is missing.")}
  if (round == TRUE & is.null(round_names)) {
    stop("Argument 'round_names' cannot be NULL if round = TRUE.")
  }

  # formats and data to start
  cat("\nPreparing data, please wait...\n\n")
  if (virtualsp_format %in% c("shp", "gpkg")) {
    if (virtualsp_format == "shp") {
      M_patt <- ".shp$"
      subs <- ".shp"
      mlist <- gsub(subs, "", list.files(path = virtualsp_folder, pattern = M_patt))
      spnames <- mlist
    } else {
      M_patt <- ".gpkg$"
      subs <- ".gpkg"
      mlist <- list.files(path = virtualsp_folder, pattern = M_patt)
      spnames <- gsub(subs, "", mlist)
    }
  } else {
    M_patt <- paste0(rformat_type(var_format), "$")
    subs <- rformat_type(var_format)
    mlist <- list.files(path = virtualsp_folder, pattern = M_patt, full.names = TRUE)
    spnames <- gsub(subs, "", list.files(path = virtualsp_folder, pattern = M_patt))
  }
  v_patt <- paste0(rformat_type(var_format), "$")

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

  bin_tabs <- lapply(1:dim(variables)[3], function(i) {
    cat("\nPreparing range values from environmental layers:\n")
    M_range <- lapply(1:length(mlist), function(j) {
      if (virtualsp_format %in% c("shp", "gpkg")) {
        if (virtualsp_format == "shp") {
          M <- rgdal::readOGR(dsn = virtualsp_folder, layer = mlist[j], verbose = FALSE)
        } else {
          M <- rgdal::readOGR(paste0(virtualsp_folder, "/", mlist[j]), spnames[j],
                              verbose = FALSE)
        }
      } else {
        M <- raster::raster(mlist[j])
      }

      mvar <- raster::mask(raster::crop(variables[[i]], M), M)
      range <- range(na.omit(mvar[]))
      cat("\t", j, "of", length(mlist), "species finished\n")

      return(range)
    })

    # overall range
    M_range <- do.call(rbind, M_range)
    overall_range <- range(c(M_range))

    if (overall_range[2] > 999) {
      overall_range <- round(overall_range / 10)
      M_range <- round(M_range / 10)
    }
    if (overall_range[2] > 9999) {
      overall_range <- round(overall_range / 100)
      M_range <- round(M_range / 100)
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
    cat("\nPreparing bin tables using ranges:\n")

    bin_table <- bin_env_null(overall_range, M_range, bin_size)
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(gsub("_", " ", spnames), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    # write table
    write.csv(bin_table,
              paste0(output_directory, "/", names(variables)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "of", dim(variables)[3], "variables processed\n")
    return(bin_table)
  })

  names(bin_tabs) <- names(variables)
  return(bin_tabs)
}
