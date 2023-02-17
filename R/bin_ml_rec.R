#' Maximum likelihood reconstruction of ancestral character states
#'
#' @param tree_data a list of two elements (phy and data) resulting from using the
#' function \code{\link[geiger]{treedata}}.
#' @param ... other arguments from \code{\link[ape]{ace}}. Arguments \code{x},
#' \code{phy}, \code{type}, and \code{method} are fixed.
#'
#' @return A table with columns representing bins, rows representing first tip
#' states and then reconstructed nodes.
#'
#' @details
#' Reconstructions are done using the function \code{\link[ape]{ace}} from the
#' \code{ape} package. The argument method is set as "ML" and the type
#' of variable is "discrete".
#'
#' @importFrom ape ace
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
#' # Maximum likelihood reconstruction
#' ml_rec <- bin_ml_rec(treeWdata)


bin_ml_rec <- function(tree_data, ...) {
  if (missing(tree_data)) {stop("Argument 'tree_data' needs to be defined.")}

  # Data from geiger::treedata
  tphy <- tree_data$phy
  ntips <- length(tphy$tip.label)
  nnode <- tphy$Nnode
  tdata <- tree_data$data

  # Matrix to fill with reconstructions
  reconMatrix <- matrix(nrow = nnode + 3, ncol = ncol(tdata))
  colnames(reconMatrix) <- colnames(tdata)
  rownames(reconMatrix) <- c(seq.int(from = 1 + ntips, to = ntips + nnode),
                             "LogLik", "Rates", "SE")

  # Reconstruct each column
  for (i in 1:ncol(tdata)) {
    # If all tips are the same, scores all the nodes for that column as the same
    if (all(tdata[, i] == tdata[1, i])) {
      reconMatrix[1:nnode, i] <- rep(tdata[1, i], nnode)
    } else{
      # Reconstruction
      temp <- ape::ace(x = tdata[, i], phy = tphy, type = "discrete",
                       method = "ML", ...)

      # Round each node to 0, 1, or ?
      alh <- temp$lik.anc
      maxlik <- round(apply(alh, 1, max), digits = 10)

      # Codes reconstructions conservatively if there is equivocation
      ancRes <- sapply(1:nnode, function(j) {
        matches <- round(alh[j, ], digits = 10) == maxlik[j]
        if (sum(matches) > 1) {return("?")} else {return(names(matches)[matches])}
      })
      reconMatrix[1:nnode,i] <- ancRes

      # Reconstruction statistics
      reconMatrix[nnode + 1, i] <- temp$loglik
      reconMatrix[nnode + 2, i] <- temp$rates
      reconMatrix[nnode + 3, i] <- temp$se
    }
  }
  whole_rec_table <- rbind(tdata, reconMatrix)

  return(whole_rec_table)
}
