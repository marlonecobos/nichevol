#' Helper function to find raster extention
#' @param format (character) any of the format types allowed for raster objects.
#' See \code{\link[raster]{writeFormats}} (e.g., "GTiff").
#' @return Raster extension according to format type.
#' @export
#' @examples
#' rformat <- rformat_type("GTiff")

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
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline layout hist plot.new points title
#' @return
#' A PDF file written in the output directory containing all resultant figures.
#' @export
#' @examples
#' # example data
#' e_data <- list(rnorm(1000, 15, 7), rnorm(800, 20, 6), rnorm(1000, 12, 3))
#' o_data <- list(sample(seq(5, 29, 0.1), 45), sample(seq(10, 33, 0.1), 40),
#'                sample(seq(1, 16, 0.1), 50))
#' y_val <- list(rep(3, length(o_data)), rep(4, length(o_data)),
#'               rep(2, length(o_data)))
#' s_names <- c("sp1", "sp2", "sp3")
#' lims <- rbind(c(3.5, 26.47), c(10.83, 29.66), c(6.92, 16.91))
#'
#' # the runing (before running, create output_directory in current directory)
#' \dontrun{
#' bins <- pdf_histograms(env_data = e_data, occ_data = o_data, y_values = y_val,
#'                        sp_names = s_names, variable_name = "Temperature",
#'                        CL_lines = 95, limits = lims, col = "geen",
#'                        output_directory = "Histogram_ranges_check")
#' }

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
#' @param bin_size (numeric) size of bins. Range of environmental values to
#' be considered when creating each character in bin tables. See details.
#' Default = 10.
#' @details
#' The argument \code{bin_size} helps to create characters that represent not
#' only one value of an environmental variable, but a range of environmental
#' conditions. For instance, if a variable of precipitation in mm is used, a
#' value of 10 for \code{bin_size} indicates that each character will represent
#' a class that correspond to 10 continuous values of precipitation (e.g., from
#' 100 to 110 mm).
#' @return
#' A character matrix (table of characters) containing bins for a given variable
#' and for all species considered. See more details in \code{\link{bin_tables}}.
#' @export
#' @examples
#' # example
#' o_range <- c(1, 25)
#' m_range <- rbind(c(5, 15), c(10, 23), c(4, 20))
#' s_range <- rbind(c(7, 15), c(12, 21), c(3, 18))
#'
#' # bin preparation
#' bins <- bin_env(overall_range = o_range, M_range = m_range,
#'                 sp_range = s_range, bin_size = 1)

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
    M_test <- seq(round(M_range[i, 1]), round(M_range[i, 2]), 1)
    sp_test <- seq(round(sp_range[i, 1]), round(sp_range[i, 2]), 1)

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
#' @param bin_size (numeric) size of bins. Range of environmental values to
#' be considered when creating each character in bin tables. See details.
#' Default = 10.
#' @details
#' The argument \code{bin_size} helps to create characters that represent not
#' only one value of an environmental variable, but a range of environmental
#' conditions. For instance, if a variable of precipitation in mm is used, a
#' value of 10 for \code{bin_size} indicates that each character will represent
#' a class that correspond to 10 continuous values of precipitation (e.g., from
#' 100 to 110 mm).
#' @export
#' @return
#' A character matrix containing bins that identify accessibility of environments
#' (environemntal values inside M; "0"), used environments (presence of species
#' records; "1"), and unknown (uncertainty about the ability of the species to be
#' in such conditions given accessibility; "?").

bin_env_null <- function(overall_range, M_range, bin_size) {
  .Deprecated("bin_env")
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
    M_test <- seq(round(M_range[i, 1]), round(M_range[i, 2]), 1)

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


#' Helper function to rename tips of trees for simulations
#' @param tree an object of class "phylo".
#' @param names (character) vector of new names. Length must be equal to number
#' of tips. They will be assigned in the order given.
#' @return Tree of class "phylo" with specified names
#' @export
#' @examples
#' # a simple tree
#' tree <- phytools::pbtree(b = 1, d = 0, n = 5, scale = TRUE,
#'                          nsim = 1, type = "continuous", set.seed(5))
#' # renaming tips
#' renamedTree <- rename_tips(tree, c("a", "b", "c", "d", "e"))

rename_tips <- function(tree, names) {
  if (missing(tree)) {stop("Argument tree needs to be defined.")}
  if (missing(names)) {stop("Argument names needs to be defined.")}
  tree$tip.label <- names
  return(tree)
}


#' Helper function to get sigma squared values for a given dataset
#' @param tree_data a list of two elements (phy and data) resulted from using
#' the function \code{\link[geiger]{treedata}}. NOTE: data must be a single
#' vector (i.e., a single column).
#' @param model model to fit to comparative data; see
#' \code{\link[geiger]{fitContinuous}}. Default = "BM".
#' @return the sigma squared value (evolutionary rate) for the data, given the
#' tree.
#' @importFrom geiger fitContinuous
#' @export
#' @examples
#' # a simple tree
#' tree <- phytools::pbtree(b = 1, d = 0, n = 5, scale = TRUE,
#'                          nsim = 1, type = "continuous", set.seed(5))
#' # simple data
#' data <- rnorm(n = length(tree$tip.label))
#' names(data) <- tree$tip.label
#' # tree with data
#' treeWdata <- geiger::treedata(tree, data)
#'
#' # Estimating sigma squared for the dataset
#' sig_sq(treeWdata)

sig_sq <- function(tree_data, model = "BM") {
  if (missing(tree_data)) {stop("Argument tree_data needs to be defined.")}
  tmp <- geiger::fitContinuous(tree_data$phy, tree_data$dat, model = model,
                               ncores = 1)
  sigmaSquared <- tmp$opt$sigsq
  return(sigmaSquared)
}


#' Helper function to calculate the median bin score for a given species
#' @param character_table data.frame containing bin scores for all species.
#' NOTE: row names must be species' names.
#' @param species_name (character) name of the species to be analyzed.
#' @param include_unknown (logical) whether or not unknown bin status should be
#' included.
#' @return Median bin value for a given species (for inferring sigma squared or
#' other comparative phylogenetic analyses requiring a single continuous variable).
#' @export
#' @examples
#' # Simulate data for single number bin labels
#' dataTable <- cbind("241" = rep("1", 5),
#'                    "242" = rep("1", 5),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#'  rownames(dataTable) <- c("GadusMorhua", "GadusMacrocephalus",
#'                           "GadusChalcogrammus", "ArctogadusGlacials",
#'                           "BoreogadusSaida")
#' # Simulate data for bin labels as strings
#' dataTableStringLabel <- cbind("241 to 244" = rep("1", 5),
#'                               "244 to 246" = c("1", "1", "0", "0", "0"),
#'                               "246 to 248" = c("1", "?", "0", "0", "0"))
#' rownames(dataTableStringLabel) <- c("GadusMorhua", "GadusMacrocephalus",
#'                                     "GadusChalcogrammus", "ArctogadusGlacials",
#'                                     "BoreogadusSaida")
#' # Use function
#' score_tip(character_table = dataTable, species_name = "GadusMorhua",
#'           include_unknown = TRUE)
#' score_tip(character_table = dataTableStringLabel, species_name = "GadusMorhua",
#'           include_unknown = FALSE)

score_tip <- function(character_table, species_name, include_unknown = FALSE) {
  if (missing(character_table)) {stop("Argument character_table needs to be defined.")}
  if (missing(species_name)) {stop("Argument species_name needs to be defined.")}
  binVals <- as.data.frame(character_table[rownames(character_table) == species_name,])
  if(include_unknown == TRUE) {
    binsWunknown <- rownames(binVals)[binVals == 1 | binVals == "?"]
    if(!grepl(x = binsWunknown[1], pattern = "to")){
      binsWunknown <- as.numeric(binsWunknown)
      score <- ((max(binsWunknown) - min(binsWunknown))/2 + min(binsWunknown))
    } else{
      bnVals <- as.numeric(unique(unlist(strsplit(binsWunknown, " to "))));
      score <- (max(bnVals) - min(bnVals))/2
    }
  }
  else{
    binsWknown <- rownames(binVals)[binVals == 1];
    if(!grepl(x = binsWknown[1], pattern = "to")){
      binsWknown <- as.numeric(binsWknown)
      score <- ((max(binsWknown) - min(binsWknown))/2 + min(binsWknown))
    } else{
      bnVals <- as.numeric(unique(unlist(strsplit(binsWknown, " to "))));
      score <- (max(bnVals) - min(bnVals))/2 + min(bnVals)
    }
  }
  return(score)
}


#' Helper function to assign bin scores to every tip in a given tree
#' @param tree_data a list of two elements (phy and data) resulting from using
#' the function \code{\link[geiger]{treedata}}.
#' @param include_unknown (logical) whether or not there are unknown tips.
#' @return a list of two elements (phy and data). Data is the median bin scored
#' as present or present + unknown.
#' @importFrom geiger treedata
#' @export
#' @examples
#' # Simulate data table
#' dataTable <- cbind("241" = rep("1", 5),
#'                    "242" = rep("1", 5),
#'                    "243" = c("1", "1", "0", "0", "0"),
#'                    "244" = c("1", "1", "0", "0", "0"),
#'                    "245" = c("1", "?", "0", "0", "0"))
#' rownames(dataTable) <- c("GadusMorhua", "GadusMacrocephalus",
#'                          "GadusChalcogrammus", "ArctogadusGlacials",
#'                          "BoreogadusSaida")
#'
#' # Simulate phylogeny
#' tree <- phytools::pbtree(b = 1, d = 0, n = 5, scale = TRUE,
#'                          nsim = 1, type = "continuous", set.seed(5))
#' tree$tip.label <- c("GadusMorhua", "GadusMacrocephalus",
#'                     "GadusChalcogrammus", "ArctogadusGlacials",
#'                     "BoreogadusSaida")
#'
#' # Unite data
#' treeWithData <- geiger::treedata(tree, dataTable)
#'
#' # Get a new tree with tips scored from median bin scores
#' score_tree(treeWithData, include_unknown = TRUE)

score_tree <- function(tree_data, include_unknown = FALSE) {
  if (missing(tree_data)) {stop("Argument tree_data needs to be defined.")}
  tCode <- unlist(lapply(tree_data$phy$tip.label, function(x) {
    score_tip(character_table = tree_data$data,
             species_name = x, include_unknown = include_unknown)
  }))
  names(tCode) <- tree_data$phy$tip.label
  twd <- geiger::treedata(tree_data$phy, tCode)
  colnames(twd$data) <- c("Median_Bin_Value")
  return(twd)
}


#' Helper function to split geographic points in 9 blocks of equal size
#' @param data matrix with longitude and latitude columns, in that order.
#' @export
#' @examples
#' # random data
#' rdata <- cbind(x = rnorm(1000, -80, 12), y = rnorm(1000, -3, 15))
#' blocks <- make_9blocks(rdata)

make_9blocks <- function(data) {
  if (missing(data)) {stop("Argument data needs to be defined.")}
  ndata <- nrow(data)
  n1 <- ceiling(ndata / 3); n2 <- n1 * 2; n3 <- ndata - n2
  xylor <- data[order(data[, 2]), ]
  grp_a <- xylor[1:n1, ]
  grp_b <- xylor[(n1 + 1):n2, ]
  grp_c <- xylor[(n2 + 1):ndata, ]
  sn1 <- ceiling(n1 / 3); sn3 <- n1 - (sn1 * 2)
  sn7 <- ceiling(n3 / 3); sn9 <- n3 - (sn7 * 2)
  all_gs <- c(rep(1, sn1), rep(2, sn1), rep(3, sn3), rep(4, sn1), rep(5, sn1),
              rep(6, sn3), rep(7, sn7), rep(8, sn7), rep(9, sn9))
  xyg <- cbind(rbind(grp_a[order(grp_a[, 1]), ], grp_b[order(grp_b[, 1]), ],
                     grp_c[order(grp_c[, 1]), ]), all_gs)
  return(xyg)
}
