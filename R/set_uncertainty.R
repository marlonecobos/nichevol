#' Set values of uncertainty towards one or both ends of the variable
#'
#' @description set_uncertainty allows to define uncertainty ("?") values
#' around values denoting presence ("1") towards one or both ends of the
#' variable in a table of binary characters.
#'
#' @param character_table a matrix of characters to represent ecological niches
#' of the species of interest. A matrix containing values "1" = presence,
#' "0" = absence, and "?" = uncertain. See \code{\link{bin_table}}.
#' @param species (character) name of the species in the table for which values
#' of uncertainty will be set.
#' @param end (character) end towards which uncertainty values ("?") will be set.
#' Options are: "high", "low", or "both".
#'
#' @details
#' Values of characters around those denoting presence ("1") are manually
#' transformed to uncertain ("?") to help producing more conservative
#' reconstructions of ancestral ecological niches. This increases uncertainty in
#' reconstructions and further niche comparisons, which reduces the events of
#' niche change that can be detected. This may be especially useful when dealing
#' with species with one or just a few known records.
#'
#' @return
#' A modified matrix of characters to represent ecological niches of the
#' species of interest.
#'
#' Potential values for characters are:
#' - "1" = the species is present in those environmental conditions.
#' - "0" = the species is not present in those environmental conditions. This is,
#' those environmental conditions inside the accessible area (M) are more extreme
#' than the ones used for the species.
#' - "?" = there is no certainty about the species presence in those environmental
#' conditions.
#'
#' @export
#'
#' @examples
#' # a character table
#' data("character_table", package = "nichevol")
#'
#' character_table[, 20:28]
#'
#' # set values of uncertainty towards the lower end of the variable for species t3
#' char_tableu <- set_uncertainty(character_table, species = "t2", end = "low")
#'
#' char_tableu[, 20:28]

set_uncertainty <- function(character_table, species, end) {
  species <- gsub("_", " ", species)
  n_col <- ncol(character_table)
  where <- which(character_table[species, ] == 1)

  if (end == "both") {
    character_table[species, character_table[species, ] != 1] <- "?"
  }
  if (end == "high") {
    if (max(where) < n_col) {
      start <- max(where) + 1
      character_table[species, start:n_col] <- "?"
    }
  }
  if (end == "low") {
    if (min(where) > 1) {
      end <- min(where) - 1
      character_table[species, 1:end] <- "?"
    }
  }

  return(character_table)
}
