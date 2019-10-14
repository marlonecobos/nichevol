#' Parsimony reconstruction of characters
#' @param whole_rec_table matrix containing all reconstructed charaters for all
#' tips and nodes. It results from using the functions \code{\link{bin_par_rec}}
#' or \code{\link{bin_ml_rec}}.
#'
#' @importFrom stats smooth
#'
#' @export

smooth_rec <- function(whole_rec_table) {
  if (missing(whole_rec_table)) {stop("Argument whole_rec_table needs to be defined.")}

  statRows <- NULL
  nrows <- nrow(whole_rec_table)

  if ("LogLik" %in% rownames(whole_rec_table)) {
    statRows <- whole_rec_table[-((nrows - 3):nrows), ]
    whole_rec_table <- whole_rec_table[1:(nrows - 3), ]
    nrows <- nrows - 3
  }

  for (k in 1:nrows){
    test <- paste0(whole_rec_table[k, ], collapse = "")
    test <- gsub("?", replacement = "u", x = test, fixed = T)

    # First, get rid of unknowns at flanks of a row that are bordered by 0s
    while (grepl(pattern = "^u+0", x = test)) {test <- gsub("u0", "00", x = test)} # At the start
    while (grepl(pattern = "0u+$", x = test)) {test <- gsub("0u", "00", x = test)} # At the end

    # Smooth estimates for unimodal reconstructions
    if (grepl(x = test, pattern = "1u+1")) {#fills in unknowns sandwiched between 1s
      while (grepl(x = test, pattern = "1u+1")) {test <- gsub("1u", "11", test)}
    }
    if (grepl(x = test, pattern = "0u+0")) {#fills in unknowns sandwiched between 0s
      while (grepl(x = test, pattern = "0u+0")) {test <- gsub("0u", "00", test)}
    }

    uSplit <- unlist(strsplit(x = test, split = "u", fixed = T))
    midString <- uSplit[nchar(uSplit) > 1] # Trims off unknowns at bin periferies
    if (grepl(pattern = "10+1", x = test)){# Forces unimodal reconstruction
      midString <- paste0(as.character(smooth(as.numeric(whole_rec_table[1, ]), )),
                          collapse = "")
    }
    test <- gsub(x = test, pattern = "[01]+", replacement = midString)
    test <- gsub("u", replacement = "?", x = test, fixed = T)
    whole_rec_table[k, ] <- unlist(strsplit(test, split = ""))
  }
  if (!is.null(statRows)){
    whole_rec_table <- rbind(whole_rec_table, statRows)
  }

  return(whole_rec_table)
}
