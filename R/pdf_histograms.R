pdf_histograms <- function(env_data, occ_data, y_values, sp_names, variable_name, CL_lines, col,
                           output_folder = "Histogram_ranges_check") {

  col1 <- rep(col, each = 2) # colors for limits in actual values

  ## frecuency of each value in M
  pdf(file = paste0(output_folder, "/Histograms_", variable_name, ".pdf"),
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
        title(main = paste0("Values and median deviations for variable ", variable_name))
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft", legend = c("", "", "  The information presented below helps to visualize the",
                                     "  distribution of values of the variable in the accessible",
                                     "  area as well as the species occurrences to facilitate the ",
                                     "  delimitation of conditions to be used in further analyses."), cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi, col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])), main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])), rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]), xlab = "Median deviation")
        points(occ_data[[x - 1]], rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })

    # next pages
    par(mfrow = c(5, 2), cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
    plots <- lapply(5:sec_length, function(x){
      ### actual values
      hist(as.numeric(names(env_data[[x]])), main = gsub("_", " ", sp_names[x]), xlab = "Variable values")
      points(as.numeric(names(occ_data[[x]])), rep(y_values[[x]][1], length(occ_data[[x]])), col = "darkgreen", pch = 1)
      abline(v = limits[x, ], col = col1)

      ### median deviation
      hist(env_data[[x]], main = gsub("_", " ", sp_names[x]), xlab = "Median deviation")
      points(occ_data[[x]], rep(y_values[[x]][2], length(occ_data[[x]])), col = "darkgreen", pch = 1)
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
        title(main = paste0("Values and median deviations for variable ", variable_name))
        legend("topleft", legend = "Description", cex = 1.1, bty = "n")
        legend("topleft", legend = c("", "", "  The information presented below helps to visualize the",
                                     "  distribution of values of the variable in the accessible",
                                     "  area as well as the species occurrences to facilitate the ",
                                     "  delimitation of conditions to be used in further analyses."), cex = 0.9, bty = "n")

        legend("topright", legend = "Symbology        ", cex = 1.1, bty = "n")
        legend("topright", legend = sym_legend, lty = lin, pch = poi, col = colss, cex = 0.9, bty = "n")

      } else {
        par(cex = 0.6, mar = (c(4.5, 4, 3.5, 1) + 0.1))
        ### actual values
        hist(as.numeric(names(env_data[[x - 1]])), main = gsub("_", " ", sp_names[x - 1]), xlab = "Variable values")
        points(as.numeric(names(occ_data[[x - 1]])), rep(y_values[[x - 1]][1], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        abline(v = limits[x - 1, ], col = col1)

        ### median deviation
        hist(env_data[[x - 1]], main = gsub("_", " ", sp_names[x - 1]), xlab = "Median deviation")
        points(occ_data[[x - 1]], rep(y_values[[x - 1]][2], length(occ_data[[x - 1]])), col = "darkgreen", pch = 1)
        limit <- floor(CL_lines * length(env_data[[x - 1]]) / 100)
        abline(v = sort(env_data[[x - 1]])[limit], col = col)
      }
    })
  }
  invisible(dev.off())
}


