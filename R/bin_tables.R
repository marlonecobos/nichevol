#' Bin tables of environmental conditions in M and for occurrences from objects
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
#' @param n_bins (numeric) number of bins to be created from the range of
#' environmental values considered when creating each character in bin tables.
#' Default = 20. See details.
#' @param bin_size (numeric) argument deprecated, use n_bins instead.
#' @param save (logical) whether or not to save the results in working directory.
#' Default = FALSE.
#' @param output_directory (character) name of the folder in which results will be
#' written.
#' @param overwrite (logical) whether or not to overwrite existing results
#' in \code{output_directory}. Default = FALSE.
#' @param verbose (logical) whether messages should be printed. Default = TRUE.
#'
#' @details
#' The percentage to be defined in \code{percentage_out} must correspond with
#' one of the confidence limits defined in \code{\link{histograms_env}}
#' (argument \code{CL_lines}). For instance, if \code{CL_lines} = 95, then
#' \code{percentage_out} can only be either 5 (keeping data inside the 95 CL) or
#' 0 (to avoid exclusion of extreme values in M).
#'
#' Excluding a certain percentage of extreme environmental values prevents the
#' algorithm from considering extremely rare environmental values in the
#' accessible area for the species (M). Being too rare, these values may have
#' never been explored by the species; therefore, including them in the process
#' of preparation of the table of characters (bin table) is risky.
#'
#' The argument \code{n_bins} helps to define how many characters (bins) will be
#' considered for the range of values in each variable. This is, a value of 20
#' determines that a range of temperature (5-25) will be split approximately
#' every 1 degree. The argument \code{bin_size} has been deprecated.
#'
#' @return
#' A list named as in \code{ranges} containing the table(s) of characters.
#' A folder named as in \code{output_directory} containing all resulting csv
#' files with the tables of characters will be created if \code{save} is set as
#' TRUE.
#'
#' Potential values for characters are:
#' - "1" = the species is present in those environmental conditions.
#' - "0" = the species is not present in those environmental conditions. This is,
#' those environmental conditions inside the accessible area (M) are more extreme
#' than the ones used for the species.
#' - "?" = there is no certainty about the species presence in those environmental
#' conditions. This happens if environmental combinations are more extreme than
#' the ones found in the accessible area (M), when environmental conditions in
#' species records are as extreme as the most extreme ones in M.
#'
#' @importFrom utils write.csv
#'
#' @export
#'
#' @usage
#' bin_tables(ranges, percentage_out = 5, n_bins = 20, bin_size, save = FALSE,
#'            output_directory, overwrite = FALSE, verbose = TRUE)
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
#' bins <- bin_tables(ranges, percentage_out = 5, n_bins = 20)
#'
#' # see arguments save and output_directory to write results in local directory

bin_tables <- function(ranges, percentage_out = 5, n_bins = 20, bin_size,
                       save = FALSE, output_directory, overwrite = FALSE,
                       verbose = TRUE) {
  # checking for potential errors
  if (missing(ranges)) {stop("Argument 'ranges' is missing.")}
  if (!missing(bin_size)) {
    warning("Argument 'bin_size' is deprecated, using 'n_bins'.")
  }
  if (save == TRUE) {
    if (missing(output_directory)) {
      stop("Argument 'output_directory' is missing.")
    } else {
      if (overwrite == FALSE & dir.exists(output_directory)) {
        stop("'output_directory' already exists, to replace it use 'overwrite' = TRUE.")
      }
      if (overwrite == TRUE & dir.exists(output_directory)) {
        unlink(x = output_directory, recursive = TRUE, force = TRUE)
      }
    }
  }

  if (verbose == TRUE) {
    message("\nPreparing bin tables using ranges:")
  }

  # directory for results
  if (save == TRUE) {dir.create(output_directory)}

  bin_tabs <- lapply(1:length(ranges), function(i) {
    # preparing ranges
    spnames <- ranges[[1]][, 1]
    cl <- paste0("M_", 100 - percentage_out, c("_lowerCL", "_upperCL"))
    sp_r <- paste0("Species_", c("lower", "upper"))

    overall_range <- range(c(ranges[[i]][, c(sp_r, cl)]))
    M_range <- ranges[[i]][, cl]
    sp_range <- ranges[[i]][, 2:3]

    # bin tables
    bin_size <- diff(overall_range) / n_bins

    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    rownames(bin_table) <- gsub("_", " ", spnames)

    # write table
    if (save == TRUE) {
      write.csv(bin_table,
                paste0(output_directory, "/", names(ranges)[i], "_bin_table.csv"),
                row.names = TRUE)
    }

    if (verbose == TRUE) {
      message(i, " of ", length(ranges), " variables processed")
    }

    return(bin_table)
  })

  names(bin_tabs) <- names(ranges)
  return(bin_tabs)
}
