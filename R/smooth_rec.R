#' Smooth character table values resulted from ancestral character state reconstructions
#'
#' @param whole_rec_table matrix containing all reconstructed characters for all
#' tips and nodes. It results from using the functions \code{\link{bin_par_rec}}
#' or \code{\link{bin_ml_rec}}.
#'
#' @return
#' The matrix of reconstructed characters with smoothed values.
#'
#' @importFrom stats smooth
#' @importFrom stringr str_extract str_replace
#' @export
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # simple matrix of data
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # ancestral reconstruction
#' parsimonyReconstruction <- bin_par_rec(treeWdata)
#'
#' # smoothing reconstructions
#' smooth_rec(parsimonyReconstruction)

smooth_rec <- function(whole_rec_table) {
  if (missing(whole_rec_table)) {stop("Argument 'whole_rec_table' needs to be defined.")}

  statRows <- NULL
  nrows <- nrow(whole_rec_table)

  if ("LogLik" %in% rownames(whole_rec_table)) {
    statRows <- whole_rec_table[(nrows - 2):nrows, ]
    whole_rec_table <- whole_rec_table[1:(nrows - 3), ]
    nrows <- nrows - 3
  }

  for (k in 1:nrows){
    test <- paste(whole_rec_table[k, ], collapse = "")
    test <- gsub("?", replacement = "u", x = test, fixed = TRUE)

    # 0s between unknowns
    if (grepl(x = test, pattern = "u0+u")) {
      while (grepl(x = test, pattern = "u0+u")) {
        pull <- stringr::str_extract(string = test, pattern = "u0+u")[1]
        pull <- gsub(unlist(strsplit(pull, split = "")), pattern = "0", replacement = "u")
        test <- stringr::str_replace(test, "u0+u", paste(pull, collapse= ""))
      }
    }

    # Unknowns between 1s
    if (grepl(x = test, pattern = "1u+1")) {#fills in unknowns sandwiched between 1s
      while (grepl(x = test, pattern = "1u+1")) {
        pull <- stringr::str_extract(string = test, pattern = "1u+1")[1]
        pull <- gsub(unlist(strsplit(pull,split = "")), pattern = "u", replacement = 1)
        test <- stringr::str_replace(test, "1u+1", paste(pull, collapse= ""))
      }
    }

    # Algorithmically smooth if there are 0s and 1s alternating, in order to yield unimodal response
    midString <- stringr::str_extract(test, "1[01]+")
    if(!is.na(midString)){
      if(nchar(midString) > 3){
        if (grepl(pattern = "10+1", x = test)){
          midString <- as.numeric(unlist(strsplit(midString, split = "")))
          midString <- paste(smooth(midString), collapse = "")
        }
        test <- stringr::str_replace(test, "1[01]+", midString)
      }
    }

    # Last check to clear any 0 or unknowns flanked by 1s
    if (grepl(x = test, pattern = "1[0u]+1")) {#fills in unknowns sandwiched between 1s
      while (grepl(x = test, pattern = "1[0u]+1")) {
        pull <- stringr::str_extract(string = test, pattern = "1[0u]+1")[1]
        pull <- gsub(unlist(strsplit(pull,split = "")), pattern = "[0u]", replacement = 1)
        test <- stringr::str_replace(test, "1[0u]+1", paste(pull, collapse= ""))
      }
    }

    # Check for 0s flanked by a 1 and an unknown
    if (grepl(x = test, pattern = "u0+1")) {
      while (grepl(x = test, pattern = "u0+1")) {
        pull <- stringr::str_extract(string = test, pattern = "u0+1")[1]
        pull <- gsub(unlist(strsplit(pull,split = "")), pattern = "u", replacement = "0")
        test <- stringr::str_replace(test, "u0+1", paste(pull, collapse= ""))
      }
    }

    if (grepl(x = test, pattern = "10+u")) {
      while (grepl(x = test, pattern = "10+u")) {
        pull <- stringr::str_extract(string = test, pattern = "10+u")[1]
        pull <- gsub(unlist(strsplit(pull,split = "")), pattern = "u", replacement = "0")
        test <- stringr::str_replace(test, "10+u", paste(pull, collapse= ""))
      }
    }

    test <- gsub("u", replacement = "?", x = test, fixed = TRUE)
    whole_rec_table[k, ] <- unlist(strsplit(test, split = ""))
  }
  if (!is.null(statRows)){
    whole_rec_table <- rbind(whole_rec_table, statRows)
  }

  return(whole_rec_table)
}
