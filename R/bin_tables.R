#' Bin tables of environemntal conditions in M and occurrences
#'
#' @description bin_tables0 helps in creating csv files with bin tables
#' of environmental conditions in M and species occurrence records. All of this
#' starting from results of previous analyses, for various species, and using
#' multiple variables.
#'
#' @param ranges list of ranges of environmental values in M and in species
#' occurrences derived from using the function \code{\link{histograms_env}}.
#' @param percentage_out (numeric) percentage of extreme environmental data in M
#' to be excluded in bin creation for further analyses. Default = 5.
#' @param bin_size (numeric) size of bins. Interval of values to be considered
#' when creating bin tables. Default = 10.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Species_E_bins".
#'
#' @importFrom utils write.csv
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all resultant csv
#' files for all variables, with bins for all species. Results will also be
#' returned as a list.

bin_tables <- function(ranges, percentage_out = 5, bin_size = 10,
                       output_directory = "Species_E_bins") {
  # checking for potential errors
  if (missing(ranges)) {stop("Argument ranges is missing.")}

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
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(as.character(ranges[[i]][, 1]), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    bin_tabs <- bin_table

    # write table
    write.csv(bin_table,
              paste0(output_directory, "/", names(ranges)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "of", length(ranges), "variables processed\n")

    return(bin_tabs)
  })

  names(bin_tabs) <- names(ranges)
  return(bin_tabs)
}
