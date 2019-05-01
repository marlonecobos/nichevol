niche_diffbars <- function(reconstructed_bins, spp_rows, present = "1", absent = "0", unknown = "0 1",
                           width = 50, height = 5, res = 300, output_dir = "Difference_bars") {

  # testing for potential errors
  if (missing(reconstructed_bins)) {
    stop("Argument reconstructed_bins is needed to perform the analyses.")
  }

  # organizing data
  nod <- as.character(reconstructed_bins[!1:nrow(reconstructed_bins) %in% spp_rows, 1])
  nodes <- gsub("\\D", "", nod)

  sps <- as.character(reconstructed_bins[1:nrow(reconstructed_bins) %in% spp_rows, 1])
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
  dir.create(output_dir) # folder

  # comparisons and plots
  comp_list <- lapply(1:nrow(all_comb), function(x){
    comp <- vector()
    for (z in 1:ncol(rec_bins)) {
      from <- rec_bins[all_comb[x, 1], z]
      to <- rec_bins[all_comb[x, 2], z]

      if (from == present & to == absent) {
        comp[z] <- "loss"
      }
      if (from == present & to == present) {
        comp[z] <- "nc"
      }
      if (from == present & to == unknown) {
        comp[z] <- "nc"
      }
      if (from == absent & to == present) {
        comp[z] <- "gain"
      }
      if (from == absent & to == absent) {
        comp[z] <- "nc"
      }
      if (from == absent & to == unknown) {
        comp[z] <- "nc"
      }
      if (from == unknown & to == absent) {
        comp[z] <- "nc"
      }
      if (from == unknown & to == present) {
        comp[z] <- "nc"
      }
      if (from == unknown & to == unknown) {
        comp[z] <- "nc"
      }
    }

    # infor for lines
    linesp <- rec_bins[all_comb[x, 1], ]

    # bar creation
    bar_name <- paste0(output_dir, "/", all_nam[x], "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    poly_lines <- sapply(1:(length(h_vertices) - 1), function(y) {
      # polygons
      if (comp[y] == "nc") {
        pcolor <- "grey90"
      } else {
        pcolor <- ifelse(comp[y] == "gain", "green1", "dodgerblue3")
      }

      xs <- c(h_vertices[y], h_vertices[y + 1], h_vertices[y + 1], h_vertices[y])

      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)

      # lines
      if (linesp[y] == unknown) {
        pcolor <- "white"
      } else {
        pcolor <- ifelse(linesp[y] == present, "grey10", "transparent")
      }

      xs <- c(h_vertices[y], h_vertices[y + 1])

      lines(x = xs, y = y_line, col = pcolor, lty = 1, lwd = 1.7)
    })
    dev.off()

  })
}
