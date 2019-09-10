#' Helper to prepare bin tables of null species

bin_env_null <- function(overall_range, M_range, bin_size) {
  # sequences
  sequence_vals <- seq(overall_range[1], overall_range[2], bin_size)

  # numeric ranges
  if (bin_size > 1) {
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
  }

  bin_tab <- lapply(1:dim(M_range)[1], function(i) {
    M_test <- seq(M_range[i, 1], M_range[i, 2], 1)

    if (bin_size > 1) {
      invar_M <- sapply(1:dim(ranges)[2], function(j) {
        if (sum(M_test %in% seq(ranges[1, j], ranges[2, j], 1)) > 0) {
          invar <- "1"
        } else {
          invar <- "?"
        }
      })
    } else {
      invar_M <- sapply(1:length(sequence_vals), function(j) {
        if (sum(M_test %in% sequence_vals[j]) > 0) {
          invar <- "1"
        } else {
          invar <- "?"
        }
      })
    }

    cat("\t", i, "of", dim(M_range)[1], "species finished\n")
    return(invar_M)
  })

  bin_tab <- do.call(rbind, bin_tab)
  if (bin_size > 1) {
    colnames(bin_tab) <- bins
  } else {
    colnames(bin_tab) <- sequence_vals
  }

  return(bin_tab)
}
