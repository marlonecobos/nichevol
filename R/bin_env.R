bin_env <- function(overall_range, M_range, sp_range, bin_size) {
  # sequences
  sequence_vals <- seq(overall_range[1], overall_range[2], bin_size)

  # numeric ranges
  ranges <- sapply(1:(length(sequence_vals) - 1), function(x){
    if (x == 1) {
      c(sequence_vals[x], sequence_vals[x + 1])
    } else {
      c(sequence_vals[x] + 1 , sequence_vals[x + 1])
    }
  })

  # character bins
  bins <- sapply(1:(length(sequence_vals) - 1), function(x){
    if (x == 1) {
      paste(sequence_vals[x], sequence_vals[x + 1], sep = " to ")
    } else {
      paste(sequence_vals[x] + 1 , sequence_vals[x + 1], sep = " to ")
    }
  })

  bin_tab <- list()

  for (i in 1:dim(M_range)[1]) {
    M_test <- seq(M_range[i, 1], M_range[i, 2], 1)
    sp_test <- seq(sp_range[i, 1], sp_range[i, 2], 1)

    invar_M <- vector()
    invar_sp <- vector()

    for (j in 1:dim(ranges)[2]) {
      if (sum(M_test %in% seq(ranges[1, j], ranges[2, j], 1)) > 0) {
        invar_M[j] <- 1
      } else {
        invar_M[j] <- 0
      }

      if (sum(sp_test %in% seq(ranges[1, j], ranges[2, j], 1)) > 0) {
        invar_sp[j] <- 100
      } else {
        invar_sp[j] <- 0
      }
    }

    invar_sum <- invar_M + invar_sp
    n <- 1:length(invar_sum)
    whereM <- range({
      if (sum(invar_sum == 1) == 0) {
        c(0, 0)
      } else {
        n[invar_sum == 1]
      }
    })
    wheresp <- range(n[invar_sum >= 100])

    places <- wheresp - whereM

    if (places[1] > 0 & whereM[1] != 0) {
      invar_sum[1:(whereM[1] - 1)] <- 1
    }

    if (places[2] <= 0) {
      invar_sum[(whereM[2] + 1):length(invar_sum)] <- 1
    }

    bin_tab[[i]] <- vector()
    for (j in 1:length(invar_sum)) {
      if(invar_sum[j] == 0) bin_tab[[i]][j] <- "?"
      if(invar_sum[j] == 1) bin_tab[[i]][j] <- "0"
      if(invar_sum[j] >= 100) bin_tab[[i]][j] <- "1"
    }

    cat("\t", i, "of", dim(M_range)[1], "species finished\n")
  }

  bin_tab <- do.call(rbind, bin_tab)
  colnames(bin_tab) <- bins

  return(bin_tab)
}
