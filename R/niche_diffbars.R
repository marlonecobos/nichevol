#' Bars for representing niche evolution
#'
#' @description nichevol_bars helps in producing bar plots that represent how
#' species niches (considering one environmental variable) have changed from
#' ancestors to decendants.
#'
#' @param reconstructed_bins matrix of reconstructed bins for nodes and species
#' derived from a process of maximum parsimony reconstruction.
#' @param species_rows (numeric) vector indicating the rows of the matrix in
#' which species bins are.
#' @param present (character) code indicating environmental bins in which the
#' species is present. Default = "1".
#' @param absent (character) code indicating environmental bins in which the
#' species is absent. Default = "0".
#' @param unknown (character) code indicating environmental bins in which the
#' species presence is unknown (uncertain). Default = "0 1".
#' @param present_col color for line representing environments where the species
#' is present. Default = "grey10".
#' @param unknown_col color for line representing environments where the species
#' presence is unknown (uncertain). Default = "white".
#' @param no_change_col color for area of the bar representing environments where
#' no change has been detected. Default = "grey90".
#' @param retraction_col color for area of the bar representing environments where
#' niche retraction has been detected. Default = "dodgerblue3".
#' @param expansion_col color for area of the bar representing environments where
#' niche expansion has been detected. Default = "green1".
#' @param width (numeric) width of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 50.
#' @param height (numeric) height of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 5.
#' @param res (numeric) nominal resolution in ppi to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 300.
#' @param output_directory (character) name of the folder in which results will
#' be written. The directory will be created as part of the process.
#' Default = "Difference_bars"
#'
#' @importFrom utils combn
#' @importFrom graphics legend lines
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all bar figures
#' produced.

nichevol_bars <- function(reconstructed_bins, species_rows, present = "1",
                          absent = "0", unknown = "0 1", present_col = "grey10",
                          unknown_col = "white", no_change_col = "grey90",
                          retraction_col = "dodgerblue3", expansion_col = "green1",
                          width = 50, height = 5, res = 300,
                          output_directory = "Difference_bars") {

  # testing for potential errors
  if (missing(reconstructed_bins)) {
    stop("Argument reconstructed_bins is needed to perform the analyses.")
  }

  # organizing data
  nod <- as.character(reconstructed_bins[!1:nrow(reconstructed_bins)
                                         %in% species_rows, 1])
  nodes <- gsub("\\D", "", nod)

  sps <- as.character(reconstructed_bins[1:nrow(reconstructed_bins)
                                         %in% species_rows, 1])
  spp <- gsub(" ", "_", sps)

  rec_bins <- reconstructed_bins[, -1]
  rec_bins <- data.frame(lapply(rec_bins, as.character), stringsAsFactors = FALSE)
  row.names(rec_bins) <- c(nod, sps)

  tpol <- ncol(rec_bins)
  wpol <- 1 / tpol

  h_vertices <- seq(0, 1, wpol)
  v_vertices <- rep(c(0, 0.05), each = 2)
  y_line <- rep(mean(c(0, 0.05)), 2)

  # combinations for comparisons
  all_comb <- t(combn(x = nod, m = 2))
  nod_sps <- cbind(rep(nod, each = length(sps)), rep(sps, length(nod)))
  all_comb <- rbind(all_comb, all_comb[, 2:1], nod_sps)

  all_nam <- t(combn(x = nodes, m = 2))
  nodes_spp <- cbind(rep(nodes, each = length(sps)), rep(spp, length(nod)))
  all_nam <- rbind(all_nam, all_nam[, 2:1], nodes_spp)
  all_nam <- paste("From_node", all_nam[, 1], "to", all_nam[, 2], sep = "_")

  # plotting data
  dir.create(output_directory) # folder

  # comparisons and plots
  comp_list <- lapply(1:nrow(all_comb), function(x){
    comp <- sapply(1:ncol(rec_bins), function(z) {
      from <- rec_bins[all_comb[x, 1], z]
      to <- rec_bins[all_comb[x, 2], z]
      if (from == present & to == absent) {comp <- "loss"}
      if (from == present & to == present) {comp <- "nc"}
      if (from == present & to == unknown) {comp <- "nc"}
      if (from == absent & to == present) {comp <- "gain"}
      if (from == absent & to == absent) {comp <- "nc"}
      if (from == absent & to == unknown) {comp <- "nc"}
      if (from == unknown & to == absent) {comp <- "nc"}
      if (from == unknown & to == present) {comp <- "nc"}
      if (from == unknown & to == unknown) {comp <- "nc"}
      return(comp)
    })

    # infor for lines
    linesp <- rec_bins[all_comb[x, 1], ]

    # bar creation
    bar_name <- paste0(output_directory, "/", all_nam[x], "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    poly_lines <- sapply(1:(length(h_vertices) - 1), function(y) {
      # polygons
      if (comp[y] == "nc") {
        pcolor <- no_change_col
      } else {
        pcolor <- ifelse(comp[y] == "gain", expansion_col, retraction_col)
      }

      xs <- c(h_vertices[y], h_vertices[y + 1], h_vertices[y + 1], h_vertices[y])
      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)

      # lines
      if (linesp[y] == unknown) {
        pcolor <- unknown_col
      } else {
        pcolor <- ifelse(linesp[y] == present, present_col, "transparent")
      }

      xs <- c(h_vertices[y], h_vertices[y + 1])
      lines(x = xs, y = y_line, col = pcolor, lty = 1, lwd = 1.7)
    })
    dev.off()

  })

  # legend
  png(filename = paste0(output_directory, "/0_Legend.png"), width = 50, height = 50,
      units = "mm", bg = "transparent", res = res)
  par(mar = rep(0, 4), cex = 1.2)
  plot(x = c(0, 0.5), y = c(0, 0.5), col = "transparent", axes = FALSE)
  legend("center",
         legend = c("Uncertain", "Present", "No change", "Retraction", "Expansion"),
         box.col = "grey94",
         lty = c(1, NA, NA, NA, NA),
         lwd = 2, col = c("transparent", NA, NA, NA, NA), bg = "grey94")

  legend("center", legend = c("                ", "", "", "", ""), bty = "n",
         pch = 22, pt.bg = c(NA, NA, no_change_col, retraction_col, expansion_col),
         pt.cex = 2.2, lty = 1, col = "transparent")

  legend("center", legend = c("                  ", "", "", "", ""), bty = "n",
         lty = c(1, 1, NA, NA, NA), lwd = 2,
         col = c(unknown_col, present_col, NA, NA, NA))
  invisible(dev.off())
}
