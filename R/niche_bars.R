niche_bars <- function(bin_table, present = "1", unknown = "?", width = 50,
                       height = 5, res = 300, output_dir = "Niche_bars") {

  # testing for potential errors
  if (missing(bin_table)) {
    stop("Argument bin_table is needed to perform the analyses.")
  }

  # organizing data
  spnames <- as.character(bin_table[, 1])
  bin_table <- bin_table[, -1]

  tpol <- ncol(bin_table)
  wpol <- 1 / tpol

  h_vertices <- seq(0, 1, wpol)
  v_vertices <- rep(c(0, 0.05), each = 2)

  dir.create(output_dir)

  barss <- sapply(1:nrow(bin_table), function(j) {
    bar_name <- paste0(output_dir, "/", gsub(" ", "_", spnames[j]), "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    polys <- sapply(1:(length(h_vertices) - 1), function(x) {
      if (as.character(bin_table[j, x]) == unknown) {
        pcolor <- "grey85"
      } else {
        pcolor <- ifelse(as.character(bin_table[j, x]) == present, "red4", "royalblue1")
      }

      xs <- c(h_vertices[x], h_vertices[x + 1], h_vertices[x + 1], h_vertices[x])

      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)
    })
    dev.off()
  })

  # legend
  png(filename = paste0(output_dir, "/0_Legend.png"), width = 50, height = 30, units = "mm",
      bg = "transparent", res = res)
  par(mar = rep(0, 4), cex = 1.2)
  plot(x = c(0, 0.5), y = c(0, 0.5), col = "transparent", axes = FALSE)
  legend("center", legend = c("Uncertain", "Present", "Not present"), box.col = "grey94",
         pch = 22, pt.bg = c("grey85", "red4", "royalblue1"), pt.cex = 2.2,
         lty = 1, col = "transparent", bg = "grey94")
  invisible(dev.off())
}



