bin_tables <- function(ranges, percentage_out = 5, bin_size = 10,
                       output_folder = "Species_E_space_bins") {

  cat("\nPreparing bin tables using ranges:\n")

  # directory for results
  dir.create(output_folder)

  bin_tabs <- list()
  for (i in 1:length(ranges)) {
    # preparing ranges
    overall_range <- range(c(ranges[[i]][, c(2:3, 6:7)]))

    cl <- paste0("M_", 100 - percentage_out, c("_lowerCL", "_upperCL"))
    M_range <- ranges[[i]][, cl]
    sp_range <- ranges[[i]][, 2:3]

    if (overall_range[2] > 999) {
      overall_range <- round(overall_range / 10)
      M_range <- round(M_range / 10)
      sp_range <- round(sp_range / 10)
    }
    if (overall_range[2] > 9999) {
      overall_range <- round(overall_range / 100)
      M_range <- round(M_range / 100)
      sp_range <- round(sp_range / 100)
    }

    # modification of range
    o_minimum <- overall_range[1]
    o_minimumc <- ifelse(o_minimum == 0, 0, floor(o_minimum / bin_size) * bin_size) - bin_size

    o_maximum <- overall_range[2]
    o_maximumc <- ifelse(o_maximum == 0, 0, ceiling(o_maximum / bin_size) * bin_size) + bin_size

    overall_range <- c(o_minimumc, o_maximumc)

    # bin tables
    bin_table <- bin_env(overall_range, M_range, sp_range, bin_size)
    bin_heads <- colnames(bin_table)
    bin_table <- data.frame(as.character(ranges[[i]][, 1]), bin_table)
    colnames(bin_table) <- c("Species", bin_heads)

    bin_tabs[[i]] <- bin_table

    # write table
    write.csv(bin_table, paste0(output_folder, "/", names(variables)[i], "_bin_table.csv"),
              row.names = FALSE)

    cat(i, "of", dim(variables)[3], "variables processed\n")
  }

  names(bin_tabs) <- names(ranges)
  return(bin_tabs)
}
