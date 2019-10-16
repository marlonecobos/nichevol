#' Parsimony reconstruction of characters
#' @param tree_data a list of two elements (phy and data) resulting from using the
#' function \code{\link[geiger]{treedata}}.
#'
#' @return A table with columns representing bins, rows representing first tip
#' states and then reconstructed nodes.
#'
#' @examples
#'
#' #Simulate data
#' tree <- phytools::pbtree(b = 1, d = 0, n = 5, scale = TRUE,
#'                          nsim = 1, type = "continuous", set.seed(5));
#' dataTable <- cbind("241" = rep("1", length(tree$tip.label)),
#'                    "242" = rep("1", length(tree$tip.label)),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"));
#' rownames(dataTable) <- tree$tip.label;
#'
#' treeWdata <- geiger::treedata(tree, dataTable);
#'
#' #Do the acutal reconstruction
#' bin_par_rec(treeWdata);
#'
#' @export

bin_par_rec <- function(tree_data) {
  if (missing(tree_data)) {stop("Argument tree_data needs to be defined.")}

  # Data fro analyses
  tphy <- tree_data$phy
  ntips <- length(tphy$tip.label)
  nnode <- tphy$Nnode
  tdata <- tree_data$data

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
      cdat <- as.integer(as.factor(tdata[,i]))
      temp <- castor::asr_max_parsimony(tree = tphy, tip_states = cdat);
      colnames(temp$ancestral_likelihoods) <- as.character(unique(tdata[,i]));

      # Round each node to 0, 1, or ? based on likelihood
      alh <- temp$ancestral_likelihoods
      maxlik <- round(apply(alh, 1, max), digits = 10)

      # Codes reconstructions conservatively if there is equivocation
      ancRes <- sapply(1:nnode, function(j) {
        matches <- round(alh[j, ], digits = 10) == maxlik[j]
        if (sum(matches) > 1) {return("?")} else {return(names(matches)[matches])}
      })
      reconMatrix[1:nnode, i] <- ancRes
    }
  }

  whole_rec_table <- rbind(tdata, reconMatrix)

  return(whole_rec_table)
}
