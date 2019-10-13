#' Parsimony reconstruction of characters
#' @param tree_data a list of two elements (phy and data) resulted from using the
#' function \code{\link[geiger]{treedata}}.
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
      temp <- castor::asr_max_parsimony(tree = tphy, tip_states = tdata[, i],
                                        Nstates = 3)
      colnames(temp$ancestral_likelihoods) <- c("1", "2", "3")

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

  # Reformat characters to follow coding for 0s, 1s, and ?s
  whole_rec_table <- gsub(whole_rec_table, pattern = "3", replacement = "?")
  whole_rec_table <- gsub(whole_rec_table, pattern = "1", replacement = "0")
  whole_rec_table <- gsub(whole_rec_table, pattern = "2", replacement = "1")

  return(whole_rec_table)
}
