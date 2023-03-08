#' Read tables of binary niche characters from directory
#'
#' @description Read one or multiple tables binary niche characters from
#' directory.
#'
#' @param file (character) name of CSV file containing a table of binary niche
#' characters.
#'
#' @rdname read_bin_table
#'
#' @return
#' A matrix if \code{read_bin_table} is used.
#'
#' A list of matrices if \code{read_bin_tables} is used.
#'
#' @export

read_bin_table <- function(file) {
  if (missing(file)) {stop("Argument 'file' must be defined.")}
  if (class(file)[1] != "character") {
    stop("'file' must be a character.")
  }

  tab <- read.csv(file[1], row.names = 1)

  cnames <- colnames(tab)

  colnames(tab) <- gsub("X", "",
                        gsub("X.", "-",
                             gsub(".to.", " to ",
                                  gsub(".to..", " to -", cnames, fixed = TRUE),
                                  fixed = TRUE),
                             fixed = TRUE))

  return(as.matrix(tab))
}


#' @param directory (character) name of directory where tables of binary niche
#' characters were written as CSV files.
#'
#' @rdname read_bin_table
#'
#' @export

read_bin_tables <- function(directory) {
  if (missing(directory)) {stop("Argument 'directory' must be defined.")}
  if (!dir.exists(directory)) {stop("'directory' does not exist.")}
  if (class(directory)[1] != "character") {
    stop("'directory' must be a character.")
  }

  files <- list.files(directory, pattern = ".csv$", full.names = TRUE)
  filenam <- list.files(directory, pattern = ".csv$")
  filenam <- gsub("_bin_table.csv", "", filenam)

  tabs <- lapply(files, read_bin_table)

  names(tabs) <- filenam

  return(tabs)
}
