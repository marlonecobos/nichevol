#' Bin tables of environemntal conditions in M and occurrences from objects
#'
#' @description bin_tables helps in creating bin tables of environmental
#' conditions in accessible areas (M) and species occurrence records
#' (i.e., table of characters). This is done using results from previous
#' analyses, and can be applied to various species and multiple variables.
#'
#' @param ranges list of ranges of environmental values in M and in species
#' occurrences derived from using the function \code{\link{histograms_env}}.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. See details. Default = 5.
#' @param bin_size (numeric) size of bins. Range of environmental values to
#' be considered when creating each character in bin tables. See details.
#' Default = 10.
#' @param save (logical) whether or not to save the results in working directory.
#' Default = FALSE.
#' @param overwrite (logical) whether or not to overwrite exitent results in
#' \code{output_directory}. Default = FALSE.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Species_E_bins".
#'
#' @details
#' The percentage to be defined in \code{percentage_out} must correspond with
#' one of the confidence limits defined in \code{\link{histograms_env}}
#' (argument \code{CL_lines}). For instance, if \code{CL_lines} = 95, then
#' \code{percentage_out} can only be either 5 (keeping data inside the 95 CL) or
#' 0 (to avoid exclussion of extreme values in M).
#'
#' Excluding a certain percentage of extreme environmental values prevents from
#' considering extremely rare environmental values in the accessible area for
#' the species (M). Being too rare, these values may have never been explored by
#' the species; therefore, including them in the process of preparation of the
#' table of characters (bin table) is risky.
#'
#' The argument \code{bin_size} helps to create characters that represent not
#' only one value of an environmental variable, but a range of environmental
#' conditions. For instance, if a variable of precipitation in mm is used, a
#' value of 10 for \code{bin_size} indicates that each character will represent
#' a class that correspond to 10 continuous values of precipitation (e.g., from
#' 100 to 110 mm).
#'
#' @return
#' A list named as in \code{ranges} containg the table(s) of characters.
#' A folder named as in \code{output_directory} containing all resultant csv
#' files with the tables of characters will be created if \code{save} is set as
#' TRUE.
#'
#' Potential values for characters are:
#'
#' "1" = the species is present in those environmental conditions.
#'
#' "0" = the species is not present in those environmental conditions. This is,
#' those environmental conditions inside the accessible area (M) are more extreme
#' than the ones used for the species.
#'
#' "?" = there is no certainty about the species presence in those environmental
#' conditions. This happens in environmental combinations more extreme than the
#' ones found in the accessible area (M), when environmental conditions in
#' species records are as extreme as the most extreme ones in M.
#'
#' @importFrom utils write.csv
#'
#' @export
#'
#' @examples
#' # simple list of ranges
#' ranges <- list(temp = data.frame(Species = c("sp1", "sp2", "sp3"),
#'                                  Species_lower = c(120, 56, 59.75),
#'                                  Species_upper = c(265, 333, 333),
#'                                  M_lower = c(93, 39, 56),
#'                                  M_upper = c(302, 333, 333),
#'                                  M_95_lowerCL = c(158, 91, 143),
#'                                  M_95_upperCL = c(292, 290, 326)),
#'                prec = data.frame(Species = c("sp1", "sp2", "sp3"),
#'                                  Species_lower = c(597, 3, 3),
#'                                  Species_upper = c(3492, 2673, 6171),
#'                                  M_lower = c(228, 3, 3),
#'                                  M_upper = c(6369, 7290, 6606),
#'                                  M_95_lowerCL = c(228, 3, 3),
#'                                  M_95_upperCL = c(3114, 2376, 2568)))
#'
#' # bin preparation
#' bins <- bin_tables(ranges, percentage_out = 5, bin_size = 10)

bin_tables <- function(ranges, percentage_out = 5, bin_size = 10, save = FALSE,
                       overwrite = FALSE, output_directory = "Species_E_bins") {
  # checking for potential errors
  if (missing(ranges)) {stop("Argument ranges is missing.")}
  if (save == TRUE) {
    if (overwrite == FALSE & dir.exists(output_directory)) {
      stop("output_directory already exists, to replace it use overwrite = TRUE.")
    }
    if (overwrite == TRUE & dir.exists(output_directory)) {
      unlink(x = output_directory, recursive = TRUE, force = TRUE)
    }
  }

  cat("\nPreparing bin tables using ranges:\n")

  # directory for results
  if (save == TRUE) {dir.create(output_directory)}

  bin_tabs <- lapply(1:length(ranges), function(i) {
    # preparing ranges
    cl <- paste0("M_", 100 - percentage_out, c("_lowerCL", "_upperCL"))
    sp_r <- paste0("Species_", c("lower", "upper"))

    overall_range <- range(c(ranges[[i]][, c(sp_r, cl)]))
    M_range <- ranges[[i]][, cl]
    sp_range <- ranges[[i]][, 2:3]

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
    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    rownames(bin_table) <- as.character(ranges[[i]][, 1])

    # write table
    if (save == TRUE) {
      write.csv(bin_table,
                paste0(output_directory, "/", names(ranges)[i], "_bin_table.csv"),
                row.names = TRUE)
    }

    cat(i, "of", length(ranges), "variables processed\n")

    return(bin_table)
  })

  names(bin_tabs) <- names(ranges)
  return(bin_tabs)
}
