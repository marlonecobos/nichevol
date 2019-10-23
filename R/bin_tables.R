#' Bin tables of environemntal conditions in M and occurrences
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
#' @param bin_size (numeric) size of bins. Interval of environmental values to
#' be considered when creating bin tables. Default = 10.
#' @param save (logical) whether or not to save the results in working directory.
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
#' @return
#' A list named as in \code{ranges} containg the table(s) of characters.
#' A folder named as in \code{output_directory} containing all resultant csv
#' files with the tables of characters will be created in \code{save} is set as
#' TRUE.
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
  dir.create(output_directory)

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
    #bin_heads <- colnames(bin_table)
    #bin_table <- data.frame(as.character(ranges[[i]][, 1]), bin_table)
    #colnames(bin_table) <- c("Species", bin_heads)
    rownames(bin_table) <- as.character(ranges[[i]][, 1])

    #bin_tabs <- bin_table

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
