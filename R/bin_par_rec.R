#' Maximum parsimony reconstruction of ancestral character states
#'
#' @param tree_data a list of two elements (phy and data) resulting from using
#' the function \code{\link[geiger]{treedata}}.
#' @param ... other arguments from \code{\link[castor]{asr_max_parsimony}}.
#' Arguments \code{tree} and \code{tip_states} are fixed.
#'
#' @return A table with columns representing bins, rows representing first tip
#' states and then reconstructed nodes.
#'
#' @details
#' Reconstructions are done using the \code{\link[castor]{asr_max_parsimony}}
#' function from the \code{castor} package.
#'
#' @importFrom ape reorder.phylo
#' @importFrom castor asr_max_parsimony
#'
#' @export
#'
#' @examples
#' # a simple tree
#' data("tree5", package = "nichevol")
#'
#' # a matrix of niche charactes (1 = present, 0 = absent, ? = unknown)
#' dataTable <- cbind("241" = rep("1", length(tree5$tip.label)),
#'                    "242" = rep("1", length(tree5$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- tree5$tip.label
#'
#' # list with two objects (tree and character table)
#' treeWdata <- geiger::treedata(tree5, dataTable)
#'
#' # Maximum parsimony reconstruction
#' par_rec <- bin_par_rec(treeWdata)


bin_par_rec <- function(tree_data, ...) {
  if (missing(tree_data)) {stop("Argument 'tree_data' needs to be defined.")}

  # Data for analyses
  tphy <- ape::reorder.phylo(tree_data$phy, order = "cladewise")
  ntips <- length(tphy$tip.label)
  nnode <- tphy$Nnode
  tdata <- tree_data$data
  tdata <- tdata[match(tphy$tip.label,row.names(tdata)),]
  tdata[tdata == 1] <- "3"
  tdata[tdata == "?"] <- "2"
  tdata[tdata == 0] <- "1"
  tdata <- apply(tdata,c(1,2),as.numeric)

  # Matrix to fill with reconstructions
  reconMatrix <- matrix(nrow = nnode, ncol = ncol(tdata))
  colnames(reconMatrix) <- colnames(tdata)
  rownames(reconMatrix) <- c(seq.int(from = 1 + ntips, to = ntips + nnode))

  # Reconstruct each column
  for (i in 1:ncol(tdata)) {
    # If all tips are the same, scores all the nodes for that column as the same
    if (all(tdata[, i] == tdata[1, i])) {
      reconMatrix[1:nnode, i] <- rep(tdata[1, i], nnode)
    } else{
      # Reconstruction
      cdat <- tdata[, i]
      temp <- castor::asr_max_parsimony(tree = tphy, tip_states = cdat, ...)
      if (length(unique(tdata[,i])) < ncol(temp$ancestral_likelihoods)){
        colnames(temp$ancestral_likelihoods) <- as.character(c(1,2,3))
      }
      else {colnames(temp$ancestral_likelihoods) <- as.character(unique(tdata[, i]))}

      # Round each node to 0, 1, or ? based on likelihood
      alh <- temp$ancestral_likelihoods
      maxlik <- round(apply(alh, 1, max), digits = 10)

      # Codes reconstructions conservatively if there is equivocation
      ancRes <- sapply(1:nnode, function(j) {
        matches <- round(alh[j, ], digits = 10) == maxlik[j]
        if (sum(matches) > 1) {return(2)} else {return(names(matches)[matches])}
      })
      reconMatrix[1:nnode, i] <- ancRes
    }
  }

  tdata <- tdata[match(tree_data$phy$tip.label,row.names(tdata)),] #Change order back to original
  whole_rec_table <- rbind(tdata, reconMatrix)

  #Converting it back to original coding
  whole_rec_table[whole_rec_table == "1"] <- "0"
  whole_rec_table[whole_rec_table == "3"] <- "1"
  whole_rec_table[whole_rec_table == "2"] <- "?"

  return(whole_rec_table)
}
