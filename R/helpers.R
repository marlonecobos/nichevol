#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}}
#' @export
#' @return Raster extension according to format type.

rformat_type <- function(format) {
  if (missing(format)) {stop("Argument format needs to be defined.")}
  if (format == "raster") {format1 <- ".grd"}
  if (format == "GTiff") {format1 <- ".tif"}
  if (format == "EHdr") {format1 <- ".bil"}
  if (format == "ascii") {format1 <- ".asc"}
  if (format == "SAGA") {format1 <- ".sdat"}
  if (format == "IDRISI") {format1 <- ".rst"}
  if (format == "CDF") {format1 <- ".nc"}
  if (format == "ENVI") {format1 <- ".envi"}
  if (format == "HFA") {format1 <- ".img"}
  return(format1)
}


#' Helper function to create PDF files with histograms
#' @param env_data list of environmental values in M for all species.
#' @param occ_data list of environmental values in occurrences for all species.
#' @param y_values list of values for the y axis to be used to represent where
#' occurrences are distributed across the environmental values in M.
#' @param sp_names (character) names of the species for which the process will
#' be performed.
#' @param variable_name (character) name of the variable to be plotted.
#' @param CL_lines (numeric) confidence limits to be plotted in the histograms.
#' @param limits numeric matrix containing the actual values for the confidence
#' limits of M.
#' @param col color for lines representing the confidence limits of M.
#' @param output_directory (character) name of the folder in which results will be
#' written. Default = "Histogram_ranges_check".
#' @export
#' @return
#' A PDF file written in the output directory containing all resultant figures.

pdf_histograms <- function(env_data, occ_data, y_values, sp_names,
                           variable_name, CL_lines, limits, col,
                           output_directory = "Histogram_ranges_check") {
  # checking for errors
  if (missing(env_data)) {stop("Argument env_data is missing.")}
  if (missing(occ_data)) {stop("Argument occ_data is missing.")}
  if (missing(y_values)) {stop("Argument y_values is missing.")}
  if (missing(sp_names)) {stop("Argument sp_names is missing.")}
  if (missing(variable_name)) {stop("Argument variable_name is missing.")}
  if (missing(CL_lines)) {stop("Argument CL_lines is missing.")}
  if (missing(limits)) {stop("Argument limits is missing.")}
  if (missing(col)) {stop("Argument col is missing.")}

  # colors for limits in actual values
  col1 <- rep(col, each = 2)

  # frecuency of each value in M
  pdf(file = paste0(output_directory, "/Histograms_", variable_name, ".pdf"),
      width = 6, height = 9)

  # legend characteristics
  sym_legend <- c("", "", "Occurrences", paste0(CL_lines, "% Confidence limits"))
  lin <- c(NA, NA, NA, rep(1, length(CL_lines)))
  poi <- c(NA, NA, 1, rep(NA, length(CL_lines)))
  colss <- c(NA, NA, "darkgreen", col)

  # plots
  if (length(env_data) > 4) {
    init_length <- 5
    sec_length <- length(env_data)

    # first page
    layout(mat = cbind(c(1, 2, 4, 6, 8), c(1, 3, 5, 7, 9)))
    plots <- lapply(1:init_length, function(x){
      if (x == 1) {
        par(cex = 0.6, mar = (c(0, 0.5, 3.5, 1.2) + 0.1))
        plot.new()
        title(main = paste0("Values and median deviations for variable ",
                            variable_name))
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft", legend = c("", "", "  The information presented below helps to visualize the",
                                     "  distribution of values of the variable in the accessible",
                                     "  area as well as the species occurrences to facilitate the ",
                                     "  delimitation of conditions to be used in further analyses."),
               cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi,
               col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])),
             main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])),
               rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])),
               col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]),
             xlab = "Median deviation")
        points(occ_data[[x - 1]],
               rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])),
               col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })

    # next pages
    par(mfrow = c(5, 2), cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
    plots <- lapply(5:sec_length, function(x){
      ### actual values
      hist(as.numeric(names(env_data[[x]])),
           main = gsub("_", " ", sp_names[x]), xlab = "Variable values")
      points(as.numeric(names(occ_data[[x]])),
             rep(y_values[[x]][1], length(occ_data[[x]])), col = "darkgreen",
             pch = 1)
      abline(v = limits[x, ], col = col1)

      ### median deviation
      hist(env_data[[x]], main = gsub("_", " ", sp_names[x]),
           xlab = "Median deviation")
      points(occ_data[[x]], rep(y_values[[x]][2], length(occ_data[[x]])),
             col = "darkgreen", pch = 1)
      limit <- floor(CL_lines * length(env_data[[x]]) / 100)
      abline(v = sort(env_data[[x]])[limit], col = col)
    })

  } else {
    init_length <- length(env_data) + 1
    # columns of layout
    c1 <- c(1, seq(2, (length(env_data) * 2), 2))
    c2 <- seq(1, ((length(env_data) * 2) + 1), 2)

    # first page
    layout(mat = cbind(c1, c2))
    plots <- lapply(1:init_length, function(x){
      if (x == 1) {
        par(cex = 0.6, mar = (c(0, 0.5, 3.5, 1.2) + 0.1))
        plot.new()
        title(main = paste0("Values and median deviations for variable ",
                            variable_name))
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft", legend = c("", "", "  The information presented below helps to visualize the",
                                     "  distribution of values of the variable in the accessible",
                                     "  area as well as the species occurrences to facilitate the ",
                                     "  delimitation of conditions to be used in further analyses."),
               cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi,
               col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])),
             main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])),
               rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])),
               col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]),
             xlab = "Median deviation")
        points(occ_data[[x - 1]],
               rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])),
               col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })
  }
  invisible(dev.off())
}


#' Helper function to prepare bin tables
#' @param overall_range (numeric) minimum and maximum values of all species and
#' Ms to be analyzed.
#' @param M_range matrix of ranges of environmental values in M for all species.
#' Columns must be minimum and maximum, and rows correspond to species.
#' @param sp_range matrix of ranges of environmental values in occurrences
#' for all species. Columns must be minimum and maximum, and rows correspond to
#' species.
#' @param bin_size (numeric) size of bins. Interval of values to be considered
#' when creating bin tables. Default = 10.
#' @export
#' @return
#' A character matrix containing bins that identify accessibility of environments
#' (environemntal values inside M; "0"), used environments (presence of species
#' records; "1"), and unknown (uncertainty about the ability of the species to be
#' in such conditions given accessibility; "?").

bin_env <- function(overall_range, M_range, sp_range, bin_size) {
  # checking for errors
  if (missing(overall_range)) {stop("Argument overall_range is missing.")}
  if (missing(M_range)) {stop("Argument M_range is missing.")}
  if (missing(sp_range)) {stop("Argument sp_range is missing.")}
  if (missing(bin_size)) {stop("Argument bin_size is missing.")}

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
    sp_test <- seq(sp_range[i, 1], sp_range[i, 2], 1)

    if (bin_size > 1) {
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
        if (sum(invar_sum == 1) == 0) {c(0, 0)} else {n[invar_sum == 1]
        }
      })
      wheresp <- range(n[invar_sum >= 100])

      places <- wheresp - whereM

      if (places[1] > 0 & whereM[1] != 0) {invar_sum[1:(whereM[1] - 1)] <- 1}

      if (places[2] <= 0) {invar_sum[(whereM[2] + 1):length(invar_sum)] <- 1}
    } else {
      invar_M <- vector()
      invar_sp <- vector()

      for (j in 1:length(sequence_vals)) {
        if (sum(M_test %in% sequence_vals[j]) > 0) {
          invar_M[j] <- 1
        } else {
          invar_M[j] <- 0
        }

        if (sum(sp_test %in% sequence_vals[j]) > 0) {
          invar_sp[j] <- 100
        } else {
          invar_sp[j] <- 0
        }
      }

      invar_sum <- invar_M + invar_sp
      n <- 1:length(invar_sum)
      whereM <- range({
        if (sum(invar_sum == 1) == 0) {c(0, 0)} else {n[invar_sum == 1]}
      })
      wheresp <- range(n[invar_sum >= 100])

      places <- wheresp - whereM

      if (places[1] > 0 & whereM[1] != 0) {invar_sum[1:(whereM[1] - 1)] <- 1}
      if (places[2] <= 0) {invar_sum[(whereM[2] + 1):length(invar_sum)] <- 1}
    }

    bin_tab <- vector()
    for (j in 1:length(invar_sum)) {
      if(invar_sum[j] == 0) bin_tab[j] <- "?"
      if(invar_sum[j] == 1) bin_tab[j] <- "0"
      if(invar_sum[j] >= 100) bin_tab[j] <- "1"
    }

    cat("\t", i, "of", dim(M_range)[1], "species finished\n")
    return(bin_tab)
  })

  bin_tab <- do.call(rbind, bin_tab)
  if (bin_size > 1) {
    colnames(bin_tab) <- bins
  } else {
    colnames(bin_tab) <- sequence_vals
  }

  return(bin_tab)
}



#' Helper function to prepare bin tables of null species
#' @param overall_range (numeric) minimum and maximum values of all species and
#' Ms to be analyzed.
#' @param M_range matrix of ranges of  environmental values in M for all species.
#' Columns must be minimum and maximum, and rows correspond to species.
#' @param bin_size (numeric) size of bins. Interval of values to be considered
#' when creating bin tables. Default = 10.
#' @export
#' @return
#' A character matrix containing bins that identify accessibility of environments
#' (environemntal values inside M; "0"), used environments (presence of species
#' records; "1"), and unknown (uncertainty about the ability of the species to be
#' in such conditions given accessibility; "?").

bin_env_null <- function(overall_range, M_range, bin_size) {
  # checking for errors
  if (missing(overall_range)) {stop("Argument overall_range is missing.")}
  if (missing(M_range)) {stop("Argument M_range is missing.")}
  if (missing(bin_size)) {stop("Argument bin_size is missing.")}

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
  if (bin_size > 1) {colnames(bin_tab) <- bins
  } else {
    colnames(bin_tab) <- sequence_vals
  }

  return(bin_tab)
}


#' Helper function to rename tips of trees
#' @param tree an object of class "phylo".
#' @param names (character) vector of new names. Length must be equal to number
#' of tips.
#' @export
rename_tips <- function(tree, names) {
  if (missing(tree)) {stop("Argument tree needs to be defined.")}
  if (missing(names)) {stop("Argument names needs to be defined.")}
  tree$tip.label <- names
  return(tree)
}


#' Helper function to get sigma squared values for a given dataset
#' @param tree_data a list of two elements (phy and data) resulted from using the
#' function \code{\link[geiger]{treedata}}.
#' @param model model to fit to comparative data; see
#' \code{\link[geiger]{fitContinuous}}. Default = "BM".
#' @export
sig_sq <- function(tree_data, model = "BM") {
  if (missing(tree_data)) {stop("Argument tree_data needs to be defined.")}
  tmp <- geiger::fitContinuous(tree_data$phy, tree_data$dat, model = model)
  sigmaSquared <- tmp$opt$sigsq
  return(sigmaSquared)
}


#' Helper function to calculate the median bin score for a given species
#' @param tree an object of class "phylo".
#' @param character_table data.frame containig characters for all species.
#' @param species_col (character) name of the column containing name of species.
#' @param species_name (character) name of the species to be analyzed.
#' @param unknown (logical) whether or not there are unknown tips.
#' @export
score_tip <- function(character_table, species_col, species_name, unknown = FALSE) {
  if (missing(character_table)) {stop("Argument character_table needs to be defined.")}
  if (missing(species_col)) {stop("Argument species_col needs to be defined.")}
  if (missing(species_name)) {stop("Argument species_name needs to be defined.")}
  binVals <- as.data.frame(character_table[character_table[, species_col] ==
                                             species_name, ])
  if(unknown == T) {
    score <- median(as.numeric(colnames(binVals)[binVals == 1 | binVals == "?"]))
  }
  else{
    score <- median(as.numeric(colnames(binVals)[binVals == 1]))
  }
  return(score)
}


#' Helper function to assign bin scores to every tip in a given tree
#' @param tree an object of class "phylo".
#' @param character_table data.frame containig characters for all species.
#' @param species_col (character) name of the column containing name of species.
#' @param unknown (logical) whether or not there are unknown tips.
#' @export
score_tree <- function(tree, character_table, species_col, unknown = FALSE) {
  if (missing(tree)) {stop("Argument species_name needs to be defined.")}
  if (missing(character_table)) {stop("Argument character_table needs to be defined.")}
  if (missing(species_col)) {stop("Argument species_col needs to be defined.")}
  tCode <- unlist(lapply(tree$tip.label, function(x) {
    binScore(character_table = character_table, species_col = species_col,
             species_name = x, unknown = unknown)
  }))
  names(tCode) <- tree$tip.label
  twd <- geiger::treedata(tree, tCode)
  return(twd)
}
