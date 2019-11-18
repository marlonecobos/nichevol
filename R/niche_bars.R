#' Bars for niche shifts among distinct taxa
#'
#' @description niche_bars produces bar figures that represent how
#' species niches have changed compared to others.
#'
#' @param bin_table matrix of environmental bins for all species derived from
#' using the functions \code{\link{bin_tables}} or \code{\link{bin_tables0}}.
#' @param present (character) code indicating environmental bins in which the
#' species is present. Default = "1".
#' @param unknown (character) code indicating environmental bins in which the
#' species presence is unknown (uncertain). Default = "?".
#' @param present_col color for line representing environments where the species
#' is present. Default = "red4".
#' @param unknown_col color for line representing environments where the species
#' presence is unknown (uncertain). Default = "grey85".
#' @param absent_col color for area of the bar representing environments where
#' no change has been detected. Default = "royalblue1".
#' @param width (numeric) width of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 50.
#' @param height (numeric) height of the device in mm to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 5.
#' @param res (numeric) nominal resolution in ppi to be passed to the
#' \code{\link[grDevices]{png}} function. Default = 300.
#' @param output_directory (character) name of the folder in which results will
#' be written. The directory will be created as part of the process.
#' Default = "Niche_bars"
#'
#' @importFrom graphics par plot polygon legend
#' @importFrom grDevices dev.off png
#'
#' @export
#'
#' @return
#' A folder named as in \code{output_directory} containing all bar figures
#' produced.

niche_bars <- function(bin_table, present = "1", unknown = "?",
                       present_col = "red4", unknown_col = "grey85",
                       absent_col = "royalblue1", width = 50, height = 5,
                       res = 300, output_directory = "Niche_bars") {

  # testing for potential errors
  if (missing(bin_table)) {
    stop("Argument 'bin_table' is needed to perform the analyses.")
  }

  # organizing data
  spnames <- as.character(bin_table[, 1])
  bin_table <- bin_table[, -1]

  tpol <- ncol(bin_table)
  wpol <- 1 / tpol

  h_vertices <- seq(0, 1, wpol)
  v_vertices <- rep(c(0, 0.05), each = 2)

  dir.create(output_directory)

  barss <- sapply(1:nrow(bin_table), function(j) {
    bar_name <- paste0(output_directory, "/", gsub(" ", "_", spnames[j]),
                       "_bar.png")

    png(filename = bar_name, width = width, height = height, units = "mm",
        bg = "transparent", res = res)
    par(mar = rep(0, 4))
    plot(x = c(0, 1), y = c(0, 0.05), col = "transparent", axes = FALSE)

    polys <- sapply(1:(length(h_vertices) - 1), function(x) {
      if (as.character(bin_table[j, x]) == unknown) {
        pcolor <- unknown_col
      } else {
        pcolor <- ifelse(as.character(bin_table[j, x]) == present,
                         present_col, absent_col)
      }

      xs <- c(h_vertices[x], h_vertices[x + 1], h_vertices[x + 1], h_vertices[x])

      polygon(x = xs, y = v_vertices, col = pcolor, border = NA)
    })
    dev.off()
  })

  # legend
  png(filename = paste0(output_directory, "/0_Legend.png"), width = 50,
      height = 30, units = "mm", bg = "transparent", res = res)
  par(mar = rep(0, 4), cex = 1.2)
  plot(x = c(0, 0.5), y = c(0, 0.5), col = "transparent", axes = FALSE)
  legend("center", legend = c("Uncertain", "Present", "Not present"),
         box.col = "grey94", pch = 22, pt.cex = 2.2, lty = 1, col = "transparent",
         pt.bg = c(unknown_col, present_col, absent_col), bg = "grey94")
  invisible(dev.off())
}



